# For AD w.r.t. to second order tensors, i.e. the deformation gradient
function Tensors._extract_gradient(v::Tuple{ForwardDiff.Dual, SVector}, t::Tensor{2, dim}) where {dim}
  P = Tensors._extract_gradient(v[1], t)
  return P, ForwardDiff.value.(v[2])
end

# For nested AD w.r.t to second order tensors
function Tensors._extract_gradient(
  v::Tuple{Tensor{2, 3, <:ForwardDiff.Dual}, SVector}, 
  t::Tensor{2, dim}
) where {dim}
  A = Tensors._extract_gradient(v[1], t)
  return A, ForwardDiff.value.(v[2])
end

# This catches the case where AD can't trace the temperature
# in hyperelastic models so we just set the entropy to zero
# since mathematically speaking dψ/dθ = 0.0 for hyperelstic
# models
function Tensors._extract_gradient(v::Tuple{T, SVector}, ::T) where T <: Number
  return -0.0, v[2]
end 

# still a work in progress
function gradient(model, props, Δt, F, θ, state_old)
  P = Tensors.gradient(z -> helmholtz_free_energy(model, props, Δt, z, θ, state_old)[1], F)
  η = -ForwardDiff.derivative(z -> helmholtz_free_energy(model, props, Δt, F, z, state_old)[1], θ)
  dstate_old = ForwardDiff.gradient(z -> helmholtz_free_energy(model, props, Δt, F, θ, z)[1], state_old)
  dstate_dstate = ForwardDiff.jacobian(z -> helmholtz_free_energy(model, props, Δt, F, θ, z)[2], state_old)
  return P, η, dstate_old, dstate_dstate
end

# fall back methods
# this way only the helmholtz_free_energy method needs to be implemented
function pk1_stress(model, props, Δt, F, θ, state_old)
  # P, state_new = Tensors.gradient(z -> helmholtz_free_energy(model, props, Δt, z, θ, state_old), F)
  # return P, state_new

  P, state_new = Tensors.gradient(z -> helmholtz_free_energy(model, props, Δt, z, θ, state_old), F)
  return P, state_new
end

function material_tangent(model, props, Δt, F, θ, state_old)
  A, state_new = Tensors.gradient(z -> pk1_stress(model, props, Δt, z, θ, state_old), F)
  return A, state_new
end

# for special cases where we just want energy and pk1_stress
# e.g. in a trust region solver
function helmholtz_free_energy_and_pk1_stress(model, props, Δt, F, θ, state_old)
  res_1, res_2 = Tensors.gradient(z -> helmholtz_free_energy(model, props, Δt, z, θ, state_old), F, :all)
  P, state_new = res_1
  ψ = Tensors._extract_value(res_2[1])
  return ψ, P, state_new
end

function pk1_stress_and_material_tangent(model, props, Δt, F, θ, state_old)
  res_1, res_2 = Tensors.gradient(z -> pk1_stress(model, props, Δt, z, θ, state_old), F, :all)
  A, state_new = res_1
  P = Tensors._extract_value(res_2[1])
  return P, A, state_new
end

function entropy(model, props, Δt, F, θ, state_old)
  η, state_new = Tensors.gradient(z -> helmholtz_free_energy(model, props, Δt, F, z, state_old), θ)
  return η, state_new
end

function mechanical_state(model, props, Δt, F, θ, state_old)
  P, ψ = Tensors.gradient(z -> helmholtz_free_energy(model, props, Δt, z, θ, state_old)[1], F, :all)
  A, state_new = Tensors.gradient(z -> pk1_stress(model, props, Δt, z, θ, state_old), F)
  return ψ, P, A, state_new
end
