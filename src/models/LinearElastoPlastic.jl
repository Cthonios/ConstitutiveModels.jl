struct LinearElastoPlastic{NP, Y, H} <: AbstractMechanicalModel{NP, 7}
  elasticity::LinearElastic
  yield_surface::Y
  isotropic_hardening::H
end

function LinearElastoPlastic(inputs::Dict{Symbol, Any})
  elasticity = LinearElastic()
  yield_surface = eval(Symbol(inputs[Symbol("yield surface")]))()
  isotropic_hardening = eval(Symbol(inputs[Symbol("isotropic hardening")]))()
  NP = num_properties(elasticity) + 
       num_properties(yield_surface) + 
       num_properties(isotropic_hardening) + 1
  return LinearElastoPlastic{NP, typeof(yield_surface), typeof(isotropic_hardening)}(
    elasticity, yield_surface, isotropic_hardening
  )
end

property_map(::LinearElastoPlastic) = Symbol[
  Symbol("density"),
  Symbol("yield surface")
]

function initialize_properties(::LinearElastoPlastic, inputs::Dict{Symbol, Any})
  elasticity = LinearElastic()
  yield_surface = eval(Symbol(inputs[Symbol("yield surface")]))()
  isotropic_hardening = eval(Symbol(inputs[Symbol("isotropic hardening")]))()
  return vcat(
    inputs[:density], 
    initialize_properties(elasticity, inputs)...,
    initialize_properties(yield_surface, inputs)...,
    initialize_properties(isotropic_hardening, inputs)...
  )
end

initialize_state(::LinearElastoPlastic) = zeros(SVector{7, Float64})

function elastic_properties(model::LinearElastoPlastic, props::V) where V <: AbstractArray{<:Number, 1}
  n_props = num_properties(model.elasticity)
  return SVector{n_props, eltype(props)}(props[1:n_props])
end

function isotropic_hardening_properties(model::LinearElastoPlastic, props::V) where V <: AbstractArray{<:Number, 1}
  # props_start = num_properties(model.elasticity) + 2
  props_start = 4
  n_yield_props = num_properties(model.yield_surface)
  n_iso_hard_props = num_properties(model.isotropic_hardening)
  n_props = n_yield_props + n_iso_hard_props
  return SVector{n_props, eltype(props)}(props[props_start:props_start + n_props - 1])
end

function unpack_state(::LinearElastoPlastic, state_old)
  ε_p_old = SymmetricTensor{2, 3, eltype(state_old), 6}(state_old.data[1:6])
  α_old   = state_old[7]
  return ε_p_old, α_old
end

function pack_state(::LinearElastoPlastic, ε_p_new, α_new)
  # ε_p_new.data
  return SVector{7, eltype(ε_p_new)}((ε_p_new.data..., α_new))
end

@inline function helmholtz_free_energy(model::LinearElastoPlastic, props, ∇u, θ, state_old, Δt)
  # extract props needed at this level
  μ              = props[3]
  elastic_props  = elastic_properties(model, props)
  iso_hard_props = isotropic_hardening_properties(model, props)

  # kinematics
  ε = symmetric(∇u)

  # unpack state variables
  ε_p_old, α_old = unpack_state(model, state_old)

  # calculate elastic trial stress
  ε_e_tr = ε - ε_p_old
  σ_e_tr = cauchy_stress(model.elasticity, elastic_props, ε_e_tr)

  # calculate hardening increment
  σ_eff = effective_stress(model.yield_surface, σ_e_tr)
  Δγ    = hardening_increment(model.isotropic_hardening, iso_hard_props, μ, σ_eff, α_old)

  # radial return
  if Δγ > 0.0
    N = dev(σ_e_tr) / norm(dev(σ_e_tr))
    ε_p_new = ε_p_old + Δγ * N
    α_new   = α_old   + sqrt(2. / 3.) * Δγ
  else
    ε_p_new = ε_p_old
    α_new   = α_old
  end

  # update stuff to calculate energy
  ε_e = ε - ε_p_new

  # energies
  ψ_e = helmholtz_free_energy(model.elasticity, elastic_props, ε_e)
  ψ_hard = energy(model.isotropic_hardening, iso_hard_props, α_new)
  ψ = ψ_e + ψ_hard

  # pack state
  state_new = pack_state(model, ε_p_new, α_new)

  return ψ, state_new
end
