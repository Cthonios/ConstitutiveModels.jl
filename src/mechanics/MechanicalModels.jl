"""
$(TYPEDEF)
"""
abstract type AbstractMechanicalModel{NP, NS} <: AbstractConstitutiveModel{NP, NS} end

include("hyperelasticity/HyperelasticModels.jl")
include("isotropic_hardening/IsotropicHardeningModels.jl")
include("yield_surfaces/YieldSurfaces.jl")

include("plasticity/PlasticityModels.jl")
