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
``
\\psi = \\frac{1}{2}\\kappa tr\\left(\\mathbf{E}\\right)^2 +
        \\mu\\mathbf{E_0}:\\mathbf{E_0}
``
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::Hencky,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    args...
)
    # unpack properties
    κ, μ = props[1], props[2]

    # kinematics
    I     = one(typeof(∇u))
    F     = ∇u + I
    E     = 0.5 * log_safe(tdot(F))
    trE   = tr(E)
    E_dev = dev(E)

    # constitutive
    ψ_vol = 0.5 * κ * trE^2
    ψ_dev = μ * dcontract(E_dev, E_dev)
    ψ     = ψ_vol + ψ_dev
    return ψ
end

# function cauchy_stress(
#     ::Hencky,
#     props, Δt,
#     ∇u, θ, Z
# )
#     # unpack properties
#     κ, μ = props[1], props[2]

#     # kinematics
#     I       = one(typeof(∇u))
#     F       = ∇u + I
#     J       = det(F)
#     # trE     = log(J)
#     E       = 0.5 * log_safe(tdot(F))
#     trE     = tr(E)
#     # E_dev   = NaNMath.pow(J, -2. / 3.) * E
#     E_dev   = dev(E)
#     # I_1_bar = tr(NaNMath.pow(J, -2. / 3.) * tdot(F))

#     # constitutive
#     # W_vol = 0.5 * κ * (0.5 * (J^2 - 1) - NaNMath.log(J))
#     # W_dev = 0.5 * μ * (I_1_bar - 3.)
#     # W_vol = 0.5 * κ * trE^2
#     # W_dev = μ * dcontract(E_dev, E_dev)
#     # ψ     = W_vol + W_dev
#     σ = (κ * trE * I + 2. * μ * E_dev) / 1.
#     Z = typeof(Z)()
#     return σ, Z
# end