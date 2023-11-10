module ConstitutiveModels

# exports
export material_tangent
export pk1_stress
export pk2_stress
export strain_energy_density

# dependencies
using ForwardDiff
using Tensors

# some top level abstract types
abstract type ConstitutiveModel end

include("MechanicalModels.jl")

end # module
