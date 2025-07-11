abstract type AbstractIsotropicHardening{NP, NS} <: AbstractConstitutiveModel{NP, NS} end

# function objective(
#     model::
# )

include("LinearIsotropicHardening.jl")
include("NoIsotropicHardening.jl")
