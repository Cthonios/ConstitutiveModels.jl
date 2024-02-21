"""
"""
abstract type AbstractMotion end
"""
"""
abstract type SimpleMotion <: AbstractMotion end
"""
Returns the deformation gradient for a given motion
"""
function deformation_gradient end
function motion_objective end

# default 
const DefaultMotionType = Tensor
deformation_gradient(motion::Type{<:SimpleMotion}, λ::T) where T <: Number = 
deformation_gradient(motion, λ)
deformation_gradient(motion::Type{<:SimpleMotion}, model::Mod, props::Props, state, λ::T) where {Mod <: MechanicalModel, Props <: AbstractArray, T <: Number} =
deformation_gradient(motion, model, props, state, λ)

include("IsochoricUniaxialStress.jl")
include("SimpleShear.jl")
include("UniaxialStrain.jl")
include("UniaxialStressDisplacementControl.jl")

export SimpleMotion
export IsochoricUniaxialStress,
       SimpleShear,
       UniaxialStrain,
       UniaxialStressDisplacementControl
       deformation_gradient
