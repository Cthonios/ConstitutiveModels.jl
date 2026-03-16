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
    return [elastic_props.κ, elastic_props.μ]
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
    ∇u::Tensor{2, 3, T, 9}, θ, Z_old, Z_new
) where T <: Number
    κ, μ = props[1], props[2]

    # kinematics
    I = one(typeof(∇u))
    F = ∇u + I
    C = tdot(F)
    λs = principal_stretchs(C)
    logλs = Tensors._map(log, λs)
    tr_logλs = logλs[1] + logλs[2] + logλs[3]
    dev_logλs = Tensors._map(x -> x - (1 / 3) * tr_logλs, logλs)

    # constitutive
    ψ_vol = 0.5κ * tr_logλs
    ψ_dev = μ * (
        dev_logλs[1] * dev_logλs[1] +
        dev_logλs[2] * dev_logλs[2] +
        dev_logλs[3] * dev_logλs[3]
    )
    return ψ_vol + ψ_dev
end

function pk1_stress(
    model::Hencky,
    props, Δt,
    ∇u::Tensor{2, 3, T, 9}, θ, Z_old, Z_new,
    ::ForwardDiffAD
) where T <: Number
    F = ∇u + one(∇u)
    S = pk2_stress(model, props, Δt, ∇u, θ, Z_old, Z_new)
    return F ⋅ S
end

function pk2_stress(
    ::Hencky,
    props, Δt,
    ∇u::Tensor{2, 3, T, 9}, θ, Z_old, Z_new
) where T <: Number
    κ, μ = props[1], props[2]

    I = one(∇u)
    F = ∇u + I

    # Right stretch tensor
    C = tdot(F)
    # U = sqrt(Symmetric(C))   # U = sqrt(C)
    logU = 0.5 * log(C)

    # Hencky logarithmic strain
    # logU = log(U)
    tr_logU = tr(logU)
    dev_logU = logU - (tr_logU / 3) * one(logU)

    # PK2 stress (second Piola-Kirchhoff)
    S_vol = κ * tr_logU * inv(C)
    S_dev = 2.0 * μ * dev_logU
    S = S_vol + S_dev
    return S
end

# function pk2_stress(
#     model::Hencky,
#     props, Δt,
#     C::SymmetricTensor{2, 3, T, 6}, θ, Z_old, Z_new
# ) where T <: Number
#     Tensors.gradient(z -> helmholtz_free_energy(model, props, Δt, z, θ, Z_old, Z_new), C)
# end

# function pk2_tangent(
#     model::Hencky,
#     props, Δt,
#     C::SymmetricTensor{2, 3, T, 6}, θ, Z_old, Z_new
# ) where T <: Number
#     Tensors.hessian(z -> helmholtz_free_energy(model, props, Δt, z, θ, Z_old, Z_new), C)
# end

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