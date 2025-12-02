struct ArrudaBoyce <: AbstractHyperelasticModel{3, 0}
end

function initialize_props(::ArrudaBoyce, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    n = inputs["n"]
    return Properties{3, eltype(elastic_props)}(
        elastic_props.κ, elastic_props.μ, sqrt(n)
    )
end

function helmholtz_free_energy(
    ::ArrudaBoyce,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    # unpack properties
    κ, μ, sqrt_n = props[1], props[2], props[3]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
    J       = det(F)
    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    I_1_bar = tr(J_m_23 * tdot(F))
    λ_chain = sqrt(I_1_bar / 3)

    β = inverse_langevin_approximation(
        λ_chain / sqrt_n, 
        TreloarApproximation() # TODO make user input
    )

    # constitutive
    ψ_vol = 0.5 * κ * (0.5 * (J^2 - 1) - log(J))
    ψ_dev = μ * sqrt_n * (β * λ_chain - sqrt_n * log(sinh(β) / β))
    ψ = ψ_vol + ψ_dev
    return ψ
end
