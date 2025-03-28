abstract type YieldSurface{NP, NS} <: AbstractMechanicalModel{NP, NS} end
yield_stress(::YieldSurface, props::V) where V <: AbstractArray = props[1]

include("VonMises.jl")
