struct FeFp{NP, Y, H} <: AbstractMechanicalModel{NP, 10}
  elasticity::Hencky
  yield_surface::Y
  isotropic_hardening::H
end

function FeFp(inputs::Dict{Symbol, Any})
  elasticity = Hencky()
  yield_surface = eval(Symbol(inputs[Symbol("yield surface")]))()
  isotropic_hardening = eval(Symbol(inputs[Symbol("isotropic hardening")]))()
  NP = num_properties(elasticity) + 
       num_properties(yield_surface) + 
       num_properties(isotropic_hardening) + 1
  return FeFp{NP, typeof(yield_surface), typeof(isotropic_hardening)}(
    elasticity, yield_surface, isotropic_hardening
  )
end

property_map(::FeFp) = Symbol[
  Symbol("density"),
  Symbol("yield surface")
]

function initialize_properties(::FeFp, inputs::Dict{Symbol, Any})
  elasticity = Hencky()
  yield_surface = eval(Symbol(inputs[Symbol("yield surface")]))()
  isotropic_hardening = eval(Symbol(inputs[Symbol("isotropic hardening")]))()
  return vcat(
    inputs[:density], 
    initialize_properties(elasticity, inputs)...,
    initialize_properties(yield_surface, inputs)...,
    initialize_properties(isotropic_hardening, inputs)...
  )
end

# initialize_state(::FeFp) = zeros(SVector{10, Float64})
initialize_state(::FeFp) = SVector{10, Float64}(one(Tensor{2, 3, Float64, 9})..., 0.0)

function elastic_properties(model::FeFp, props::V) where V <: AbstractArray{<:Number, 1}
  n_props = num_properties(model.elasticity)
  return SVector{n_props, eltype(props)}(props[1:n_props])
end

function isotropic_hardening_properties(model::FeFp, props::V) where V <: AbstractArray{<:Number, 1}
  # props_start = num_properties(model.elasticity) + 2
  props_start = 4
  n_yield_props = num_properties(model.yield_surface)
  n_iso_hard_props = num_properties(model.isotropic_hardening)
  n_props = n_yield_props + n_iso_hard_props
  return SVector{n_props, eltype(props)}(props[props_start:props_start + n_props - 1])
end

function unpack_state(::FeFp, state_old)
  # ε_p_old = SymmetricTensor{2, 3, eltype(state_old), 6}(state_old.data[1:6])
  F_p_old = Tensor{2, 3, eltype(state_old), 9}((state_old.data[1:9]))
  α_old   = state_old[10]
  return F_p_old, α_old
end

function pack_state(::FeFp, F_p_new, α_new)
  # ε_p_new.data
  return SVector{10, eltype(F_p_new)}((F_p_new.data..., α_new))
end

@inline function helmholtz_free_energy(model::FeFp, props, ∇u, θ, state_old, Δt)
  # extract props needed at this level
  μ              = props[3]
  elastic_props  = elastic_properties(model, props)
  iso_hard_props = isotropic_hardening_properties(model, props)

  # kinematics
  F = deformation_gradient(∇u)

  # unpack state variables
  F_p_old, α_old = unpack_state(model, state_old)
  F_e_tr = dot(F, inv(F_p_old))
  # C_e_tr = right_cauchy_green(F_e_tr)
  # E_e_tr = hencky_strain(C_e_tr)

  # calculate elastic trial stress
  # really the mandel stress below
  M_e_tr = cauchy_stress(model.elasticity, elastic_props, F_e_tr - one(F_e_tr))
  # display(M_e_tr)
  # display(norm(M_e_tr))
  # calculate hardening increment
  M_eff = effective_stress(model.yield_surface, M_e_tr)
  Δγ    = hardening_increment(model.isotropic_hardening, iso_hard_props, μ, M_eff, α_old)

  # radial return
  if Δγ > 0.0
    @show "here"
    N = dev(M_e_tr) / norm(dev(M_e_tr)) |> symmetric
    # ε_p_new = ε_p_old + Δγ * N
    display(N)
    display(expm(Δγ * N))
    F_p_new = dot(exp(Δγ * N), F_p_old)
    α_new   = α_old   + sqrt(2. / 3.) * Δγ

    display(F_p_new)
  else
    # @show "hur"
    F_p_new = F_p_old
    α_new   = α_old
  end
  # display(F)
  # update stuff to calculate energy
  # ε_e = ε - ε_p_new
  # @show typeof(F)
  # @show typeof(F_p_new)
  @show size(F)
  @show size(F_p_new)  
  F_e_new = dot(F, inv(F_p_new))
  # C_e_new = right_cauchy_green(F_e_new)
  # E_e_new = hencky_strain(C_e_new)

  # energies
  # display(F_e_new)
  # display(one(F_e_new))
  ψ_e = helmholtz_free_energy(model.elasticity, elastic_props, F_e_new - one(F_e_new))
  ψ_hard = energy(model.isotropic_hardening, iso_hard_props, α_new)
  ψ = ψ_e + ψ_hard
  
  # pack state
  state_new = pack_state(model, F_p_new, α_new)

  return ψ, state_new
end
