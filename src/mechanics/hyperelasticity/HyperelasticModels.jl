"""
$(TYPEDEF)
"""
abstract type AbstractHyperelasticModel{NP, NS} <: AbstractMechanicalModel{NP, NS} end

include("Gent.jl")
include("Hencky.jl")
include("LinearElastic.jl")
include("NeoHookean.jl")
