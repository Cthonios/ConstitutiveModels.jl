"""
$(TYPEDEF)
"""
struct LinearElastic <: AbstractHyperelasticModel{2, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::LinearElastic, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    return [elastic_props.λ, elastic_props.μ]
end

"""
``\\psi = \\frac{1}{2}\\lambda tr\\left(\\varepsilon\\right)^2
        + \\mu \\varepsilon:\\varepsilon``
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::LinearElastic,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    # unpack properties
    λ, μ = props[1], props[2]

    # kinematics
    ε = linear_strain(∇u)

    # constitutive
    ψ = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)

    return ψ
end

function cauchy_stress(
    model::LinearElastic,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    # kinematics
    ε = linear_strain(∇u)
    # constitutive
    return cauchy_stress(model, props, Δt, ε, θ, Z_old, Z_new)
end

function pk1_stress(
    model::LinearElastic, 
    props, Δt, 
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    F = ∇u + one(∇u)
    J = det(F)
    σ = cauchy_stress(model, props, Δt, ∇u, θ, Z_old, Z_new)
    P = J * dot(σ, inv(F)')
    return P
end

function cauchy_stress(
    ::LinearElastic,
    props, Δt,
    ε::SymmetricTensor{2, 3, T, 6}, θ, Z_old, Z_new
) where T <: Number
    # unpack properties
    λ, μ = props[1], props[2]
    # constitutive
    I = one(ε)
    σ = λ * tr(ε) * I + 2. * μ * ε
    return σ
end
