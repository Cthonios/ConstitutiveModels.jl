"""
$(TYPEDEF)
"""
abstract type AbstractHyperelasticModel{NP, NS} <: AbstractMechanicalModel{NP, NS} end

include("ArrudaBoyce.jl")
include("Gent.jl")
include("Hencky.jl")
include("LinearElastic.jl")
include("NeoHookean.jl")
