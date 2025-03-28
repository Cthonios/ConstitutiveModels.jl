module ConstitutiveModels

using DiffResults
using DocStringExtensions
using FiniteDiff
using ForwardDiff
using MuladdMacro
using StaticArrays
using Tensors

include("AbstractTypes.jl")
# include("ContinuumTensors.jl")
include("kinematics/Eigen.jl")
include("kinematics/Kinematics.jl")
include("kinematics/SimpleMotions.jl")

# AD stuff, should we make this optional?
include("AD.jl")

include("ScalarSolver.jl")

# kernels
include("hyperelasticity/Hyperelasticity.jl")
include("isotropic_hardening/IsotropicingHardenings.jl")
include("yield_surfaces/YieldSurfaces.jl")

# models
include("models/FeFp.jl")
include("models/Hyperelastic.jl")
include("models/LinearElastoPlastic.jl")
# include("models/SimpleFeFv.jl")

# Main exposed types
export ConstitutiveModel
export initialize_properties
export initialize_state

# helpers for development
export MaterialState

export IsochoricUniaxialStress
export SimpleShear
export UniaxialStrain
export UniaxialStressDisplacementControl

# models
export FeFp
export Hyperelastic
export LinearElastoPlastic

# methods
export cauchy_stress
export deformation_gradient
export displacement_gradient
export helmholtz_free_energy
export material_tangent
export pk1_stress

end # module
