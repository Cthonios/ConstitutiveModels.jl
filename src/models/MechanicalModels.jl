abstract type MechanicalModel{NP, NS} <: ConstitutiveModel{NP, NS} end
abstract type HyperelasticModel{NP, NS} <: MechanicalModel{NP, NS} end
abstract type PlasticityModel{NP, NS} <: MechanicalModel{NP, NS} end

# include("hyperelasticity/Hencky.jl") # still needs more work on log strain
include("hyperelasticity/LinearElastic.jl")
include("hyperelasticity/NeoHookean.jl")
include("isotropic_hardening/IsotropicingHardenings.jl")
include("yield_surfaces/YieldSurfaces.jl")

include("plasticity/LinearElastoPlasticity.jl")

function MechanicalModel(
  model::Type{M}, inputs::D
) where {M <: MechanicalModel, D <: Dict}
  return model(inputs)
end

export MechanicalModel

export LinearElastic,
       NeoHookean