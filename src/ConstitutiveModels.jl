module ConstitutiveModels

# exports
export ConstitutiveModel
export initialize_state
export number_of_properties
export number_of_state_variables

export MechanicalModel
export HyperelasticModel
export LinearElastic
export NeoHookean

export cauchy_stress
export energy
export pk1_stress
export pk2_stress

# motions
export IsochoricUniaxialStress
export SimpleMotion
export SimpleShear
export UniaxialStrain
export UniaxialStressDisplacementControl
export deformation_gradient

# dependencies
using DocStringExtensions
using NonlinearSolve
using StaticArrays
using Tensors

# for docs
@template (FUNCTIONS, METHODS, MACROS) = 
"""
$(TYPEDSIGNATURES)
$(DOCSTRING)
$(METHODLIST)
"""

@template (TYPES) = 
"""
$(TYPEDFIELDS)
$(DOCSTRING)
"""

# some top level abstract types
abstract type ConstitutiveModel{NProps, NStateVars} end
function initialize_state end
number_of_properties(::ConstitutiveModel{NProps, NStateVars}) where {NProps, NStateVars} = NProps
number_of_state_variables(::ConstitutiveModel{NProps, NStateVars}) where {NProps, NStateVars} = NStateVars
Base.size(::ConstitutiveModel{NProps, NStateVars}) where {NProps, NStateVars} = (NProps, NStateVars)

include("MechanicalModels.jl")
include("Motions.jl")

end # module
