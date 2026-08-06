"""
$(TYPEDEF)
"""
struct LinearElasticity <: AbstractHyperelasticity
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::LinearElasticity, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    return [elastic_props.λ, elastic_props.μ]
end

num_properties(::LinearElasticity) = 2

function property_names(::LinearElasticity)
    return ["Lamé's first constant", "shear modulus"]
end

"""
``\\psi = \\frac{1}{2}\\lambda tr\\left(\\varepsilon\\right)^2
        + \\mu \\varepsilon:\\varepsilon``
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(model::LinearElasticity, props, ∇u::Tensor{2, 3, T, 9}, θ) where T
    ε = linear_strain(∇u)
    return helmholtz_free_energy(model, props, ε, θ)
end

function helmholtz_free_energy(::LinearElasticity, props, ε::SymmetricTensor{2, 3, T, 6}, θ) where T
    λ, μ = props[1], props[2]
    ψ = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)
    return ψ
end

function cauchy_stress(model::LinearElasticity, props, ∇u::Tensor{2, 3, T, 9}, θ) where T
    ε = linear_strain(∇u)
    return cauchy_stress(model, props, ε, θ)
end

function cauchy_stress(::LinearElasticity, props, ε::SymmetricTensor{2, 3, T, 6}, θ) where T <: Number
    λ, μ = props[1], props[2]
    I = one(ε)
    σ = λ * tr(ε) * I + 2. * μ * ε
    return σ
end

function pk1_stress(model::LinearElasticity, props, ∇u, θ)
    F = ∇u + one(∇u)
    J = det(F)
    σ = cauchy_stress(model, props, ∇u, θ)
    P = J * dot(σ, inv(F)')
    return P
end

p_wave_modulus(::LinearElasticity, props) = props[1] + 2 * props[2]

function spatial_tangent(model::LinearElasticity, props, ∇u::Tensor{2, 3, T, 9}, θ) where T
    ε = linear_strain(∇u)
    return spatial_tangent(model, props, ε, θ)
end

function spatial_tangent(::LinearElasticity, props, ::SymmetricTensor{2, 3, T, 6}, θ) where T
    λ, μ = props[1], props[2]
    δ(i, j) = i == j ? 1.0 : 0.0
    f = (i, j, k, l) -> λ * δ(i, j) * δ(k, l) + 
        μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
    return SymmetricTensor{4, 3, T, 36}(f)
end
