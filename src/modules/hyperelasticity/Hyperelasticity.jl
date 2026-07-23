"""
$(TYPEDEF)
"""
abstract type AbstractHyperelasticity <: AbstractConstitutiveModule end
num_state_variables(::AbstractHyperelasticity) = 0

# expected minimum interface to be defined by new models
function helmholtz_free_energy(::AbstractHyperelasticity, props, ∇u, θ)
    throw(UnImplementedMethodError("Need to implement helmholtz_free_energy for hyperelasticity model"))
end

function p_wave_modulus(::AbstractHyperelasticity, props)
    throw(UnImplementedMethodError("Need to implement p_wave_modulus for hyperelasticity model"))
end

# AD default
function pk1_stress(model::AbstractHyperelasticity, props, ∇u, θ)
    return Tensors.gradient(z -> helmholtz_free_energy(model, props, z, θ), ∇u)
end

# AD default
function material_tangent(::AbstractHyperelasticity, props, ∇u, θ)
    return Tensors.gradient(z -> pk1_stress(model, props, z, θ), ∇u)
end

# derived interface
function cauchy_stress(model::AbstractHyperelasticity, props, ∇u, θ)
    F = ∇u + one(∇u)
    J = det(F)
    P = pk1_stress(model, props, ∇u, θ)
    σ = (1 / J) * dot(P, transpose(F))
    return σ
end

# todo other tangent methods

include("ArrudaBoyce.jl")
include("Gent.jl")
include("Hencky.jl")
include("LinearElasticity.jl")
include("MooneyRivlin.jl")
include("NeoHookean.jl")
include("SaintVenantKirchhoff.jl")
include("SethHill.jl")
