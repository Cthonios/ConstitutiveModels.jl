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

# thermo stuff below as fallbacks

function pk1_stress(model::T, props, ∇u, θ, state_old, Δt) where T <: ConstitutiveModel
  P, state = Tensors.gradient(z -> helmholtz_free_energy(model, props, z, θ, state_old, Δt), ∇u)
  return P, state
end

@inline function material_tangent(model, props, ∇u, θ, state_old, Δt)
  A, state = Tensors.gradient(y -> 
    Tensors.gradient(z -> helmholtz_free_energy(model, props, z, θ, state_old, Δt), y),
    ∇u
  )
  return A, state
end
