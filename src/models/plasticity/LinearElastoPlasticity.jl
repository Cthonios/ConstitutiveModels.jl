struct LinearElastoPlasticity{NP, Yield, Hard} <: PlasticityModel{NP, 7}
  elastic_model::LinearElastic
  yield_surface::Yield
  isotropic_hardening_model::Hard
end

function LinearElastoPlasticity(inputs::Dict{String, Any})
  yield_surf_type = eval(Meta.parse(inputs["yield surface"]))
  iso_hard_type = eval(Meta.parse(inputs["isotropic hardening model"]))

  elastic_model = LinearElastic(inputs)
  yield_surf_model = yield_surf_type(inputs)
  iso_hard_model = iso_hard_type(inputs)

  n_props = num_properties(elastic_model) + 
            num_properties(yield_surf_model) + 
            num_properties(iso_hard_model)
  return LinearElastoPlasticity{
    n_props, typeof(yield_surf_model), typeof(iso_hard_model)
  }(elastic_model, yield_surf_model, iso_hard_model)
end

function LinearElastoPlasticity(inputs::Dict{Symbol, Any})
  yield_surf_type = eval(Meta.parse(inputs[Symbol("yield surface")]))
  iso_hard_type = eval(Meta.parse(inputs[Symbol("isotropic hardening model")]))

  elastic_model = LinearElastic(inputs)
  yield_surf_model = yield_surf_type(inputs)
  iso_hard_model = iso_hard_type(inputs)

  n_props = num_properties(elastic_model) + 
            num_properties(yield_surf_model) + 
            num_properties(iso_hard_model)
  return LinearElastoPlasticity{
    n_props, typeof(yield_surf_model), typeof(iso_hard_model)
  }(elastic_model, yield_surf_model, iso_hard_model)
end

# function LinearElastoPlasticity_old(inputs::D, type::Type{ArrType}) where {D <: Dict{String, Any}, ArrType}
#   yield_surf_type = eval(Meta.parse(inputs["yield surface"]))
#   iso_hard_type = eval(Meta.parse(inputs["isotropic hardening model"]))

#   elastic_model, elastic_props, elastic_state = LinearElastic(inputs, type)
#   yield_surf_model, yield_surf_props, yield_surf_state = yield_surf_type(inputs, type)
#   iso_hard_model, iso_hard_props, iso_hard_state = iso_hard_type(inputs, type)

#   props = vcat(elastic_props, yield_surf_props, iso_hard_props)
#   state = vcat(elastic_state, yield_surf_state, iso_hard_state)

#   return LinearElastoPlasticity{
#     num_properties(elastic_model) + num_properties(yield_surf_model) + num_properties(iso_hard_model),
#     typeof(yield_surf_model), typeof(iso_hard_model)
#   }(
#     elastic_model, yield_surf_model, iso_hard_model
#   ), props, state
# end

function initialize_props(model::LinearElastoPlasticity, inputs::Dict{Symbol, Any})
  elastic_props = initialize_props(model.elastic_model, inputs)
  yield_surf_props = initialize_props(model.yield_surface, inputs)
  iso_hard_props = initialize_props(model.isotropic_hardening_model, inputs)
  return vcat(elastic_props, yield_surf_props, iso_hard_props)
end 

function elastic_properties(model::LinearElastoPlasticity, props::V) where V <: AbstractArray{<:Number, 1}
  n_props = num_properties(model.elastic_model)
  return SVector{n_props, eltype(props)}(props.data[1:n_props])
end

function isotropic_hardening_properties(model::LinearElastoPlasticity, props::V) where V <: AbstractArray{<:Number, 1}
  props_start = num_properties(model.elastic_model) + 1
  n_yield_props = num_properties(model.yield_surface)
  n_iso_hard_props = num_properties(model.isotropic_hardening_model)
  n_props = n_yield_props + n_iso_hard_props
  return SVector{n_props, eltype(props)}(props.data[props_start:props_start + n_props - 1])
end

# function helmholtz_free_energy(model::LinearElastoPlasticity, props, ε::SymmetricTensor{2, 3, T1, 6}, state_old::SVector{7, T2}) where {T1, T2}

#   # extract props needed at this level
#   μ              = props[2]
#   elastic_props  = elastic_properties(model, props)
#   iso_hard_props = isotropic_hardening_properties(model, props)

#   # # kinematics
#   # I = one(SymmetricTensor{2, 3, eltype(F), 6})
#   # ∇u = F - I
#   # ε = symmetric(∇u)

#   # unpack state variables
#   ε_p_old = @views fromvoigt(SymmetricTensor{2, 3, eltype(ε), 6}, state_old[1:6])
#   α_old   = state_old[7]

#   # calculate elastic trial stress
#   ε_e_tr = ε - ε_p_old
#   σ_e_tr = cauchy_stress(model.elastic_model, elastic_props, ε_e_tr)

#   # calculate hardening increment
#   σ_eff = effective_stress(model.yield_surface, σ_e_tr)
#   Δγ    = hardening_increment(model.isotropic_hardening_model, iso_hard_props, μ, σ_eff, α_old)

#   # radial return
#   if Δγ > 0.0
#     N = dev(σ_e_tr) / norm(dev(σ_e_tr))
#     ε_p_new = ε_p_old + Δγ * N
#     α_new   = α_old   + sqrt(2. / 3.) * Δγ
#   else
#     ε_p_new = ε_p_old
#     α_new   = α_old
#   end

#   # update stuff to calculate energy
#   ε_e = ε - ε_p_new
#   ψ = helmholtz_free_energy(model.elastic_model, props, ε_e) + 
#       energy(model.isotropic_hardening_model, props, α_new)
  
#   return ψ
# end

function unpack_state(::LinearElastoPlasticity, state_old)
  ε_p_old = SymmetricTensor{2, 3, eltype(state_old), 6}(state_old.data[1:6])
  α_old   = state_old[7]
  return ε_p_old, α_old
end

function pack_state(::LinearElastoPlasticity, ε_p_new, α_new)
  # ε_p_new.data
  return SVector{7, eltype(ε_p_new)}((ε_p_new.data..., α_new))
end

function helmholtz_free_energy(
  model::LinearElastoPlasticity, 
  props, Δt, F, θ, state_old
)
  # extract props needed at this level
  μ              = props[2]
  elastic_props  = elastic_properties(model, props)
  iso_hard_props = isotropic_hardening_properties(model, props)

  # kinematics
  I = one(SymmetricTensor{2, 3, eltype(F), 6})
  ∇u = F - I
  ε = symmetric(∇u)

  # unpack state variables
  ε_p_old, α_old = unpack_state(model, state_old)

  # calculate elastic trial stress
  ε_e_tr = ε - ε_p_old
  σ_e_tr = cauchy_stress(model.elastic_model, elastic_props, ε_e_tr)

  # calculate hardening increment
  σ_eff = effective_stress(model.yield_surface, σ_e_tr)
  Δγ    = hardening_increment(model.isotropic_hardening_model, iso_hard_props, μ, σ_eff, α_old)

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
  ψ_e = helmholtz_free_energy(model.elastic_model, elastic_props, ε_e)
  ψ_hard = energy(model.isotropic_hardening_model, iso_hard_props, α_new)
  ψ = ψ_e + ψ_hard

  # pack state
  state_new = pack_state(model, ε_e, α_new)

  return ψ, state_new
end

# function cauchy_stress(model::LinearElastoPlasticity, props, F::Tensor{2, 3, T, 9}, state_old::SVector{7, T}) where T

#   # extract props needed at this level
#   μ              = props[2]
#   elastic_props  = elastic_properties(model, props)
#   iso_hard_props = isotropic_hardening_properties(model, props)

#   # kinematics
#   I = one(SymmetricTensor{2, 3, eltype(F), 6})
#   ∇u = F - I
#   ε = symmetric(∇u)

#   # unpack state variables
#   ε_p_old = @views fromvoigt(SymmetricTensor{2, 3, eltype(ε), 6}, state_old[1:6])
#   α_old   = state_old[7]

#   # calculate elastic trial stress
#   ε_e_tr = ε - ε_p_old
#   σ_e_tr = cauchy_stress(model.elastic_model, elastic_props, ε_e_tr)

#   # calculate hardening increment
#   σ_eff = effective_stress(model.yield_surface, σ_e_tr)
#   Δγ    = hardening_increment(model.isotropic_hardening_model, iso_hard_props, μ, σ_eff, α_old)

#   # radial return
#   if Δγ > 0.0
#     N = dev(σ_e_tr) / norm(dev(σ_e_tr))
#     σ = σ_e_tr - 2. * μ * Δγ * N
#     ε_p_new = ε_p_old + Δγ * N
#     α_new   = α_old   + sqrt(2. / 3.) * Δγ
#   else
#     σ = σ_e_tr
#     ε_p_new = ε_p_old
#     α_new   = α_old
#   end

#   # pack new state variables
#   state_new = SVector{7, eltype(ε_p_new)}(vcat(tovoigt(SVector, ε_p_new), α_new))

#   return σ_e_tr, props, state_new
# end

# function cauchy_stress(model::LinearElastoPlasticity, props, ε::SymmetricTensor{2, 3, T, 6}, state_old::SVector{7, T}) where T

#   # extract props needed at this level
#   μ              = props[2]
#   elastic_props  = elastic_properties(model, props)
#   iso_hard_props = isotropic_hardening_properties(model, props)

#   # kinematics
#   # I = one(SymmetricTensor{2, 3, eltype(F), 6})
#   # ∇u = F - I
#   # ε = symmetric(∇u)

#   # unpack state variables
#   ε_p_old = @views fromvoigt(SymmetricTensor{2, 3, eltype(ε), 6}, state_old[1:6])
#   α_old   = state_old[7]

#   # calculate elastic trial stress
#   ε_e_tr = ε - ε_p_old
#   σ_e_tr = cauchy_stress(model.elastic_model, elastic_props, ε_e_tr)

#   # calculate hardening increment
#   σ_eff = effective_stress(model.yield_surface, σ_e_tr)
#   Δγ    = hardening_increment(model.isotropic_hardening_model, iso_hard_props, μ, σ_eff, α_old)

#   # radial return
#   if Δγ > 0.0
#     N = dev(σ_e_tr) / norm(dev(σ_e_tr))
#     σ = σ_e_tr - 2. * μ * Δγ * N
#     ε_p_new = ε_p_old + Δγ * N
#     α_new   = α_old   + sqrt(2. / 3.) * Δγ
#   else
#     σ = σ_e_tr
#     ε_p_new = ε_p_old
#     α_new   = α_old
#   end

#   # pack new state variables
#   state_new = SVector{7, eltype(ε_p_new)}(vcat(tovoigt(SVector, ε_p_new), α_new))

#   return σ_e_tr, props, state_new
# end

# function pk1_stress(model::LinearElastoPlasticity, props, F, state_old)
#   return cauchy_stress(model, props, F, state_old)
# end
