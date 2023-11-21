# super type
abstract type MechanicalModel{NProps, NStateVars} <: ConstitutiveModel{NProps, NStateVars} end

# special modules that can be used on their own or recycled below

# hyperelasticity
abstract type HyperelasticModel{NProps, NStateVars} <: MechanicalModel{NProps, NStateVars} end

initialize_state(model::Mod, ::Type{Vector{<:T}}) where {Mod <: HyperelasticModel, T <: Number} = 
zeros(T, number_of_state_variables(model))
# zeros(MVector{number_of_state_variables(model), type})
function initialize_state(model::Mod, ::Type{T}) where {Mod <: HyperelasticModel, T <: Union{SVector, MVector}}
  if T <: SVector
    @info "Note SVectors can't be used with Enzyme"
  end
  state = zeros(T)
  @assert length(state) == number_of_state_variables(model)
  return state
end

# modles to include
# include("hyperelastic_models/LinearElastic.jl")
include("hyperelastic_models/NeoHookean.jl")

# abstract type PlasticityModel{NProps, NStateVars} <: MechanicalModel{NProps, NStateVars} end
# include("plasticity_models/YieldSurfaces.jl")
# include("plasticity_models/HardeningModels.jl")
# include("plasticity_models/LinearElastoPlasticity.jl")
