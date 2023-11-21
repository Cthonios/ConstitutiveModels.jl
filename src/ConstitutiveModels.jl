module ConstitutiveModels

# exports
export ConstitutiveModel
export initialize_state
export number_of_properties
export number_of_state_variables

export MechanicalModel
export HyperelasticModel
export NeoHookean

export cauchy_stress
export energy
export pk1_stress
# export pk1_stress_and_property_adjoints
export pk2_stress

# motions
export IsochoricUniaxialStress
export SimpleMotion
export SimpleShear
export UniaxialStrain
export UniaxialStressDisplacementControl
export deformation_gradient

# states
export MaterialState
export MechanicalState
# export update_deformation_gradients
export update!

# dependencies
using DocStringExtensions
using NonlinearSolve
using StaticArrays
using Tensors
using Unitful

# some top level abstract types
abstract type ConstitutiveModel{NProps, NStateVars} end
function initialize_state end
number_of_properties(::ConstitutiveModel{NProps, NStateVars}) where {NProps, NStateVars} = NProps
number_of_state_variables(::ConstitutiveModel{NProps, NStateVars}) where {NProps, NStateVars} = NStateVars
Base.size(::ConstitutiveModel{NProps, NStateVars}) where {NProps, NStateVars} = (NProps, NStateVars)

include("Uitls.jl")
include("MechanicalModels.jl")
include("Motions.jl")

include("MaterialStates.jl")

end # module
