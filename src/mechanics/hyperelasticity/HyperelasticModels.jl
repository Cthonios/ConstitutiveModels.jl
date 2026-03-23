"""
$(TYPEDEF)
All methods have the following type signature
(::AbstractHyperelastic, props, F::DeformationGradient, θ)

Models that call this module are responsible for first calculating
the deformation gradient from the displacement gradient
"""
abstract type AbstractHyperelastic{NP, NS} <: AbstractMechanicalModule{NP, NS} end

# defaults
function cauchy_stress(model::AbstractHyperelastic, props, F, θ)
    J = det(F)
    P = pk1_stress(model, props, F, θ)
    return (1 / J) * symmetric(P ⋅ transpose(F))
end
function material_tangent(model::AbstractHyperelastic, props, F, θ)
    return material_tangent(model, props, F, θ, ForwardDiffAD())
end
function material_tangent(model::AbstractHyperelastic, props, F, θ, ::ForwardDiffAD)
    return Tensors.gradient(z -> pk1_stress(model, props, z, θ), F)
end
function pk1_stress(model::AbstractHyperelastic, props, F, θ)
    return pk1_stress(model, props, F, θ, ForwardDiffAD())
end
function pk1_stress(model::AbstractHyperelastic, props, F, θ, ::ForwardDiffAD)
    return Tensors.gradient(z -> strain_energy_density(model, props, z, θ), F)
end
# needs to be implemented, no fallback
function strain_energy_density end

# meta hyperelastic model
struct Hyperelastic{
    NP, NS,
    StrainEnergyDensity <: AbstractHyperelastic{NP, NS}
} <: AbstractMechanicalModel{NP, NS}
    strain_energy_density::StrainEnergyDensity
end

function initialize_model(::Type{<:Hyperelastic}, inputs::Dict{String})
    strain_energy_density = eval(Symbol(inputs["strain energy density"]))
    return Hyperelastic(strain_energy_density())
end

function initialize_props(model::Hyperelastic, inputs::Dict{String})
    return initialize_props(model.strain_energy_density, inputs)
end

function cauchy_stress(
    model::Hyperelastic,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    return cauchy_stress(model.strain_energy_density, props, F, θ)
end

function helmholtz_free_energy(
    model::Hyperelastic,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    return strain_energy_density(model.strain_energy_density, props, F, θ)
end

function material_tangent(
    model::Hyperelastic,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    return material_tangent(model.strain_energy_density, props, F, θ)
end

function pk1_stress(
    model::Hyperelastic,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    return pk1_stress(model.strain_energy_density, props, F, θ)
end

include("ArrudaBoyce.jl")
include("Gent.jl")
include("Hencky.jl")
include("LinearElastic.jl")
include("MooneyRivlin.jl")
include("NeoHookean.jl")
include("SaintVenantKirchhoff.jl")
include("SethHill.jl")
