"""
$(TYPEDEF)
"""
struct Hencky <: AbstractHyperelasticModel{2, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::Hencky, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    return Properties{2, eltype(elastic_props)}(
        elastic_props.κ, elastic_props.μ
    )
end

"""
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::Hencky,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    # unpack properties
    κ, μ = props[1], props[2]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
    J       = det(F)
    # trE     = log(J)
    E       = 0.5 * log(tdot(F))
    trE     = tr(E)
    # devE    = dev(E)
    E_dev   = dev(E)
    # E_dev   = NaNMath.pow(J, -2. / 3.) * E
    # I_1_bar = tr(NaNMath.pow(J, -2. / 3.) * tdot(F))

    # constitutive
    # W_vol = 0.5 * κ * (0.5 * (J^2 - 1) - NaNMath.log(J))
    # W_dev = 0.5 * μ * (I_1_bar - 3.)
    W_vol = 0.5 * κ * trE^2
    W_dev = μ * dcontract(E_dev, E_dev)
    ψ     = W_vol + W_dev
    Z     = typeof(Z)()
    return ψ, Z
end

"""
$(TYPEDSIGNATURES)
"""
function cauchy_stress(
    ::Hencky,
    props, Δt,
    ∇u, θ, Z
)
    # unpack properties
    κ, μ = props[1], props[2]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
    J       = det(F)
    # trE     = log(J)
    E       = 0.5 * log_safe(tdot(F))
    trE     = tr(E)
    # E_dev   = NaNMath.pow(J, -2. / 3.) * E
    E_dev   = dev(E)
    # I_1_bar = tr(NaNMath.pow(J, -2. / 3.) * tdot(F))

    # constitutive
    # W_vol = 0.5 * κ * (0.5 * (J^2 - 1) - NaNMath.log(J))
    # W_dev = 0.5 * μ * (I_1_bar - 3.)
    # W_vol = 0.5 * κ * trE^2
    # W_dev = μ * dcontract(E_dev, E_dev)
    # ψ     = W_vol + W_dev
    σ = (κ * trE * I + 2. * μ * E_dev) / 1.
    Z = typeof(Z)()
    return σ, Z
end