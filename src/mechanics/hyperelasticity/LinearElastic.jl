"""
$(TYPEDEF)
"""
struct LinearElastic <: AbstractHyperelastic{2, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::LinearElastic, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    return [elastic_props.λ, elastic_props.μ]
end

function strain_energy_density(model::LinearElastic, props, F::Tensor{2, 3, T, 9}, θ) where T <: Number
    ∇u = F - one(F)
    ε = linear_strain(∇u)
    return strain_energy_density(model, props, ε, θ)
end

"""
``\\psi = \\frac{1}{2}\\lambda tr\\left(\\varepsilon\\right)^2
        + \\mu \\varepsilon:\\varepsilon``
$(TYPEDSIGNATURES)
"""
function strain_energy_density(::LinearElastic, props, ε::SymmetricTensor{2, 3, T, 6}, θ) where T <: Number
    # unpack properties
    λ, μ = props[1], props[2]

    # constitutive
    ψ = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)

    return ψ
end

function cauchy_stress(model::LinearElastic, props, F, θ)
    # kinematics
    ∇u = F - one(F)
    ε = linear_strain(∇u)
    # constitutive
    return cauchy_stress(model, props, ε, θ)
end

function cauchy_stress(::LinearElastic, props, ε::SymmetricTensor{2, 3, T, 6}, θ) where T <: Number
    # unpack properties
    λ, μ = props[1], props[2]
    # constitutive
    I = one(ε)
    σ = λ * tr(ε) * I + 2. * μ * ε
    return σ
end

function pk1_stress(model::LinearElastic, props, F, θ)
    ∇u = F - one(F)
    ε = linear_strain(∇u)
    J = det(F)
    σ = cauchy_stress(model, props, ε, θ)
    P = J * σ ⋅ inv(F)'
    return P
end
