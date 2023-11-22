module ConstitutiveModels

# exports
export ConstitutiveModel
export initialize_state
export number_of_properties
export number_of_state_variables

export MechanicalModel
export HyperelasticModel
export LinearElastic
export LinearElastoPlasticity
export NeoHookean

export cauchy_stress
export energy
export pk1_stress

# motions
export IsochoricUniaxialStress
export SimpleMotion
export SimpleShear
export UniaxialStrain
export UniaxialStressDisplacementControl
export deformation_gradient

# dependencies
using DocStringExtensions
using LinearSolve
using NonlinearSolve
using StaticArrays
using Tensors
import Tensors: tdot, dott

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

# some math helpers to be consistent with Tensors.jl
Tensors.tdot(F::Matrix{T}) where T <: Number = tr(F' * F)
Tensors.tdot(F::MMatrix{3, 3, T, 9}) where T <: Number = tr(F' * F)
Tensors.tdot(F::SMatrix{3, 3, T, 9}) where T <: Number = tr(F' * F)
Tensors.dott(F::Matrix{T}) where T <: Number = tr(F * F')
Tensors.dott(F::MMatrix{3, 3, T, 9}) where T <: Number = tr(F * F')
Tensors.dott(F::SMatrix{3, 3, T, 9}) where T <: Number = tr(F * F')

include("MechanicalModels.jl")
include("Motions.jl")

end # module
