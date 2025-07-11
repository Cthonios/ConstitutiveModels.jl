module ConstitutiveModels

# exports
export cauchy_stress,
       entropy,
       heat_capacity,
       heat_flux,
       helmholtz_free_energy,
       initialize_props,
       initialize_state,
       material_hessian,
       material_tangent,
       pk1_stress
export Gent,
       Hencky,
       LinearElastic, 
       NeoHookean
export LinearIsotropicHardening,
       NoIsotropicHardening
export TrescaYieldSurface,
       VonMisesYieldSurface
export LinearElastoPlasticity
export LinearThermoElastic
export deformation_gradient,
       displacement_gradient,
       simulate_material_point
export BiaxialStrain,
       IsochoricUniaxialStress,
       PureShearStrain,
       SimpleShear,
       UniaxialStrain,
       UniaxialStressDisplacementControl

# dependencies
using DocStringExtensions
using ForwardDiff
using MuladdMacro
using NaNMath
using Roots
using StaticArrays
using Tensors

include("utils/AD.jl")
include("utils/Eigen.jl")
include("utils/ElasticConstants.jl")
include("utils/solvers/NewtonSolver.jl")
include("Interface.jl")

include("mechanics/MechanicalModels.jl")
include("thermal/ThermalModels.jl")
include("thermomechanical/ThermomechanicalModels.jl")

# 
include("CommonMethods.jl")

# some testing utils
include("utils/SimpleMotions.jl")

end