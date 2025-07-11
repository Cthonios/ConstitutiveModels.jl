abstract type AbstractYieldSurface{
    NP, 
    NS,
    IH <: AbstractIsotropicHardening
} <: AbstractConstitutiveModel{NP, NS} end

function effective_stress(::AbstractYieldSurface, ::SymmetricTensor{2, 3, T, 6}) where T <: Number
    @assert "Needs to be implemented"
end

include("TrescaYieldSurface.jl")
include("VonMisesYieldSurface.jl")
