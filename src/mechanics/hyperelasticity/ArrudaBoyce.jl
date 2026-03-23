function _dinverse_langevin_approximation(x, type)
    return ForwardDiff.derivative(z -> inverse_langevin_approximation(z, type), x)
end

struct ArrudaBoyce <: AbstractHyperelastic{3, 0}
end

function initialize_props(::ArrudaBoyce, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    n = inputs["n"]
    return [elastic_props.κ, elastic_props.μ, sqrt(n)]
end

@inline function _treloar_inverse_langevin(x)
    return x * (3 - x^2) / (1 - x^2)
end

@inline function _dtreloar_inverse_langevin(x)
    return ForwardDiff.derivative(_treloar_inverse_langevin, x)
end

function strain_energy_density(::ArrudaBoyce, props, F, θ)
    κ, μ, sqrt_n = props[1], props[2], props[3]
    J       = det(F)
    Jm13    = one(eltype(∇u)) / cbrt(J)
    Jm23    = Jm13 * Jm13
    I1bar   = tr(Jm23 * tdot(F))
    λ_chain = sqrt(I1bar / 3)
    β       = _treloar_inverse_langevin(λ_chain / sqrt_n)
    ψ_vol   = 0.5κ * (0.5 * (J^2 - 1) - log(J))
    ψ_dev   = μ * sqrt_n * (β * λ_chain - sqrt_n * log(sinh(β) / β))
    return ψ_vol + ψ_dev
end

function pk1_stress(::ArrudaBoyce, props, F, θ)
    κ, μ, sqrt_n = props[1], props[2], props[3]
    J       = det(F)
    Jm13    = one(eltype(∇u)) / cbrt(J)
    Jm23    = Jm13 * Jm13
    I1bar   = tr(Jm23 * tdot(F))
    λ_chain = sqrt(I1bar / 3)
    β       = _treloar_inverse_langevin(λ_chain / sqrt_n)
    FinvT   = inv(F)'
    P_vol   = (κ / 2) * (J^2 - 1) * FinvT
    coeff   = μ * sqrt_n * β * Jm23 / (3 * λ_chain)
    P_dev   = coeff * (F - (I1bar / 3) * FinvT)
    return P_vol + P_dev
end

function material_tangent(::ArrudaBoyce, props, F, θ)
    κ, μ, sqrt_n = props[1], props[2], props[3]

    J       = det(F)
    J2      = J * J
    Jm23    = cbrt(J)^(-2)
    C       = tdot(F)
    IC      = inv(C)
    I2      = one(IC)
    trC     = tr(C)
    I1bar   = Jm23 * trC
    λ_chain = sqrt(I1bar / 3)
    β       = _treloar_inverse_langevin(λ_chain / sqrt_n)
    dβ_dλ   = _dtreloar_inverse_langevin(λ_chain / sqrt_n) / sqrt_n

    α     = μ * sqrt_n * β * Jm23 / (3 * λ_chain)
    V     = I2 - (I1bar / 3) * IC

    # PK2
    S = (κ / 2) * (J2 - 1) * IC + α * V

    # standard building blocks
    ICxIC  = IC ⊗ IC
    ICodIC = 0.5 * (otimesu(IC, IC) + otimesl(IC, IC))

    # ── ℂ_vol: identical to NeoHookean ───────────────────────────────────────
    ℂ_vol = κ * (J2 * ICxIC - (J2 - 1) * ICodIC)

    # ── ℂ_dev = 2∂(αV)/∂C ────────────────────────────────────────────────────
    dα_dλ = μ * sqrt_n * Jm23 / 3 * (dβ_dλ * λ_chain - β) / λ_chain^2
    dα2   = (dα_dλ * Jm23 / (3 * λ_chain)) * V - (2 * α / 3) * IC
    dV2   = -(2/3) * (I2 ⊗ IC) + (2 * I1bar / 3) * ICodIC

    ℂ_dev = dα2 ⊗ V + α * dV2

    ℂ = ℂ_vol + ℂ_dev
    return _convect_tangent(ℂ, S, F)
end

# function helmholtz_free_energy(
#     ::ArrudaBoyce,
#     props, Δt,
#     ∇u, θ, Z_old, Z_new
# )
#     # unpack properties
#     κ, μ, sqrt_n = props[1], props[2], props[3]

#     # kinematics
#     I       = one(typeof(∇u))
#     F       = ∇u + I
#     J       = det(F)
#     J_m_13  = 1. / cbrt(J)
#     J_m_23  = J_m_13 * J_m_13
#     I_1_bar = tr(J_m_23 * tdot(F))
#     λ_chain = sqrt(I_1_bar / 3)

#     β = inverse_langevin_approximation(
#         λ_chain / sqrt_n, 
#         TreloarApproximation() # TODO make user input
#     )

#     # constitutive
#     ψ_vol = 0.5 * κ * (0.5 * (J^2 - 1) - log(J))
#     ψ_dev = μ * sqrt_n * (β * λ_chain - sqrt_n * log(sinh(β) / β))
#     ψ = ψ_vol + ψ_dev
#     return ψ
# end

# function pk1_stress(
#     ::ArrudaBoyce,
#     props, Δt,
#     ∇u, θ, Z_old, Z_new
# )
#     κ, μ, sqrt_n = props[1], props[2], props[3]

#     F       = ∇u + one(typeof(∇u))
#     J       = det(F)
#     Jm13    = one(eltype(∇u)) / cbrt(J)
#     Jm23    = Jm13 * Jm13
#     I1bar   = tr(Jm23 * tdot(F))
#     λ_chain = sqrt(I1bar / 3)
#     β       = inverse_langevin_approximation(λ_chain / sqrt_n, TreloarApproximation())
#     FinvT   = inv(F)'

#     P_vol = (κ / 2) * (J^2 - 1) * FinvT
#     coeff = μ * sqrt_n * β * Jm23 / (3 * λ_chain)
#     P_dev = coeff * (F - (I1bar / 3) * FinvT)

#     return P_vol + P_dev
# end
# function pk1_stress(
#     ::ArrudaBoyce,
#     props, Δt,
#     ∇u, θ, Z_old, Z_new
# )
#     @show "Hur"
#     κ, μ, sqrt_n = props[1], props[2], props[3]

#     F       = ∇u + one(typeof(∇u))
#     J       = det(F)
#     Jm13    = one(eltype(∇u)) / cbrt(J)
#     Jm23    = Jm13 * Jm13
#     I1bar   = tr(Jm23 * tdot(F))
#     λ_chain = sqrt(I1bar / 3)
#     β       = inverse_langevin_approximation(λ_chain / sqrt_n, TreloarApproximation())
#     FinvT   = inv(F)'

#     P_vol = (κ / 2) * (J^2 - 1) * FinvT
#     coeff = μ * sqrt_n * β * Jm23 / (3 * λ_chain)
#     P_dev = coeff * (F - (I1bar / 3) * FinvT)

#     return P_vol + P_dev
# end

p_wave_modulus(::ArrudaBoyce, props) = props[1] + 4 * props[2] / 3
