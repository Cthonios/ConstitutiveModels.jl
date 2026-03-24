module ConstitutiveModels

# exports
export pow
export cauchy_stress,
       entropy,
       heat_capacity,
       heat_flux,
       helmholtz_free_energy,
       initialize_model,
       initialize_props,
       initialize_state,
       state_variable_names,
       p_wave_modulus,
       material_hessian,
       material_tangent,
       pk1_stress
# hyperelastic
export ArrudaBoyce,
       Gent,
       Hencky,
       LinearElastic,
       MooneyRivlin,
       NeoHookean,
       SaintVenantKirchhoff,
       SethHill
# isotropic hardening
export LinearIsotropicHardening,
       NoIsotropicHardening
# yield surfaces
export VonMisesYieldSurface
# meta models
export Hyperelastic,
       LinearElastoPlasticity,
       FiniteDefJ2Plasticity
export LinearThermoElastic
# motions
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
using InverseLangevinApproximations
using Roots
using StaticArrays
using Tensors

include("utils/Eigen.jl")
include("utils/ElasticConstants.jl")
include("utils/Kinematics.jl")
include("utils/MaterialSymmetry.jl")
include("utils/TensorUtils.jl")
include("utils/solvers/NewtonSolver.jl")

include("AbstractTypes.jl")
include("CommonMethods.jl")
include("SimpleMotions.jl")

# modules below

# conduction
include("modules/conduction/FouriersLaw.jl")
# hyperelasticity
include("modules/hyperelasticity/ArrudaBoyce.jl")
include("modules/hyperelasticity/Gent.jl")
include("modules/hyperelasticity/Hencky.jl")
include("modules/hyperelasticity/LinearElastic.jl")
include("modules/hyperelasticity/MooneyRivlin.jl")
include("modules/hyperelasticity/NeoHookean.jl")
include("modules/hyperelasticity/SaintVenantKirchhoff.jl")
include("modules/hyperelasticity/SethHill.jl")
# isotropic hardening
include("modules/isotropic_hardening/LinearIsotropicHardening.jl")
include("modules/isotropic_hardening/NoIsotropicHardening.jl")
# viscosity
include("modules/viscosity/Quadratic.jl")
# yield surfaces
include("modules/yield_surfaces/VonMisesYieldSurface.jl")

# actual models below
include("models/FiniteDefJ2Plasticity.jl")
include("models/Hyperelastic.jl")
include("models/LinearElastoPlasticity.jl")
include("models/LinearThermoElastic.jl")

# testing utils
# include("utils/SimpleMotions.jl")

end