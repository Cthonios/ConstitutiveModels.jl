# super type
abstract type MechanicalModel{NProps, NStateVars} <: ConstitutiveModel{NProps, NStateVars} end

function initialize_properties(model::Mod, props_in::V; type::Type = MVector) where {Mod <: MechanicalModel, V <: AbstractArray{<:Number, 1}}
  @assert length(props_in) == number_of_properties(model)
  if type <: MVector
    props = MVector{number_of_properties(model), eltype(V)}(props_in)
  elseif type <: SVector
    props = SVector{number_of_properties(model), eltype(V)}(props_in)
  elseif type <: Vector
    props = props_in
  else
    error("Bad type in initialize_properties for $(supertype(typeof(model))) model type.")
  end
  return props
end

# special modules that can be used on their own or recycled below

# hyperelasticity
abstract type HyperelasticModel{NProps, NStateVars} <: MechanicalModel{NProps, NStateVars} end



function initialize_state(model::Mod; type::Type = MVector) where Mod <: HyperelasticModel
  if type <: Vector
    state = zeros(T, 0)
  elseif type <: MVector
    state = zeros(MVector{number_of_state_variables(model), Float64})
  elseif type <: SVector
    @warn "Note these can't be used with in place methods"
    @warn "Note SVectors can't be used with Enzyme"
    # TODO make type here parameteric
    state = zeros(SVector{number_of_state_variables(model), Float64})
  else
    @warn "This likely isn't supported"
  end
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
