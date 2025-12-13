abstract type AbstractIsotropicHardening{NP, NS} <: AbstractConstitutiveModel{NP, NS} end

include("LinearIsotropicHardening.jl")
include("NoIsotropicHardening.jl")
