abstract type YieldSurface{NP, NS} <: ConstitutiveModel{NP, NS} end
yield_stress(::YieldSurface, props::V) where V <: AbstractArray = props[1]

include("J2YieldSurface.jl")
