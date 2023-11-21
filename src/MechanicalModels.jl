abstract type MechanicalModel{NProps, NStateVars} <: ConstitutiveModel{NProps, NStateVars} end

abstract type HyperelasticModel{NProps, NStateVars} <: MechanicalModel{NProps, NStateVars} end

function Base.zeros(type::Type, mod::Mod) where Mod <: HyperelasticModel
  n_props      = number_of_properties(mod)
  n_state_vars = number_of_state_variables(mod)
  return SVector{n_props, type}, SVector{n_state_vars, type}
end

initialize_state(model::Mod, type::Type = Float64) where Mod <: HyperelasticModel = 
zeros(SVector{number_of_state_variables(model), type})

update_state(::Mod, state::SVector{NStateVars, <:Number}) where {Mod <: HyperelasticModel, NStateVars} =
state

# modles to include
include("hyperelastic_models/LinearElastic.jl")
include("hyperelastic_models/NeoHookean.jl")

abstract type PlasticityModel{NProps, NStateVars} <: MechanicalModel{NProps, NStateVars} end
include("plasticity_models/YieldSurfaces.jl")
include("plasticity_models/HardeningModels.jl")
include("plasticity_models/LinearElastoPlasticity.jl")
