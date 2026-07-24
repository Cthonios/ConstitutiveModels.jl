module ConstitutiveModels

# methods
export as_tensor,
       cauchy_stress,
       cauchy_stress_temperature_modulus,
       deformation_gradient,
       density,
       displacement_gradient,
       dissipation,
       entropy,
       heat_capacity,
       heat_flux,
       helmholtz_free_energy,
       initialize_props,
       initialize_state,
       linear_strain,
       material_tangent,
       num_properties,
       num_state_variables,
       pk1_stress,
       pk2_stress,
       pow,
       p_wave_modulus,
       simulate_material_point,
       spatial_tangent,
       state_variable_names,
       velocity_gradient
# heat conduction models
export FouriersLaw
# hyperelasticity models
export ArrudaBoyce,
       Gent,
       Hencky,
       LinearElasticity,
       MooneyRivlin,
       NeoHookean,
       SaintVenantKirchhoff,
       SethHill
# isotropic hardening models
export LinearIsotropicHardening,
       NoIsotropicHardening,
       Voce
# material symmetry types
export Isotropy
# objective stress rates
export TruesdellCauchyStressRate
# simple motions for testing
export BiaxialStrain,
       IsochoricUniaxialStress,
       PureShearStrain,
       SimpleShear,
       UniaxialStrain,
       UniaxialStressDisplacementControl
# thermal expansion
export LinearThermalExpansion
# yield surfaces
export VonMises

# actual models
export FiniteDefJ2Plasticity
export Hyperelastic
export Hypoelastic
export LinearElastic
export LinearElastoplastic
export LinearThermoelastic
export LinearViscoelastic
export Thermal

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
include("utils/StateVariables.jl")
include("utils/TensorUtils.jl")

include("Interface.jl")

# modules
include("modules/heat_capacity/HeatCapacity.jl")
include("modules/heat_conduction/HeatConduction.jl")
include("modules/hyperelasticity/Hyperelasticity.jl")
include("modules/objective_stress_rates/ObjectiveStressRates.jl")
include("modules/plasticity/Plasticity.jl")
include("modules/thermal_expansion/ThermalExpansion.jl")

# models
include("models/FiniteDefJ2Plasticity.jl")
include("models/Hyperelastic.jl")
include("models/Hypoelastic.jl")
include("models/LinearElastic.jl")
include("models/LinearElastoplastic.jl")
include("models/LinearThermoelastic.jl")
include("models/LinearViscoelastic.jl")
include("models/Thermal.jl")

# some testing utils
# TODO should we move this elsewhere?
include("utils/SimpleMotions.jl")

end