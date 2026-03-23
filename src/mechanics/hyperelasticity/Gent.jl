"""
$(TYPEDEF)
"""
struct Gent <: AbstractHyperelastic{3, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::Gent, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    Jm = inputs["Jm"]
    return [elastic_props.κ, elastic_props.μ, Jm]
end

"""
``\\psi = \\frac{1}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
        - \\frac{1}{2}\\mu J_m\\ln\\left(1 - \\frac{\\bar{I}_1 - 3}{Jm}\\right)``
$(TYPEDSIGNATURES)
"""
function strain_energy_density(::Gent, props, F, θ)
    # unpack properties
    κ, μ, Jm = props[1], props[2], props[3]

    # kinematics
    J       = det(F)
    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    I_1_bar = tr(J_m_23 * tdot(F))

    # constitutive
    ψ_vol = 0.5 * κ * (0.5 * (J * J - 1) - log(J))
    ψ_dev = -μ * Jm / 2. * log(1. - (I_1_bar - 3.) / Jm)
    ψ     = ψ_vol + ψ_dev
    return ψ
end

function pk1_stress(::Gent, props, F, θ)
    # unpack properties
    κ, μ, Jm = props[1], props[2], props[3]

    # kinematics
    J       = det(F)
    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    I_1     = tr(tdot(F))
    I_1_bar = J_m_23 * I_1
    F_inv_T = inv(F)'

    # constitutive
    P = 0.5 * κ * (J * J - 1.) * F_inv_T + 
        μ * J_m_23 * Jm * (F - (1. / 3.) * I_1 * F_inv_T) / (Jm - I_1_bar + 3)
    return P
end

function material_tangent(::Gent, props, F, θ)
    κ, μ, Jm = props[1], props[2], props[3]

    J    = det(F)
    J2   = J * J
    Jm23 = cbrt(J)^(-2)
    C    = tdot(F)
    IC   = inv(C)
    I2   = one(IC)
    trC  = tr(C)
    β    = Jm - Jm23 * trC + 3
    γ    = μ * Jm * Jm23
    γβ   = γ / β
    V    = I2 - (trC / 3) * IC

    # PK2
    S = (κ / 2) * (J2 - 1) * IC + γβ * V

    # standard building blocks
    ICxIC  = IC ⊗ IC
    ICodIC = 0.5 * (otimesu(IC, IC) + otimesl(IC, IC))

    # same as NeoHookean
    ℂ_vol = κ * (J2 * ICxIC - (J2 - 1) * ICodIC)

    # ℂ_dev = 2∂S_dev/∂C
    # 2∂(γ/β)/∂C = 2(γ/β)*(-1/3 C⁻¹ + J^{-2/3}/β * V)
    dγβ2 = 2 * γβ * (-(1/3) * IC + (Jm23 / β) * V)

    # 2∂V/∂C = -2/3*(I⊗C⁻¹) + 2*trC/3 * C⁻¹⊙C⁻¹
    dV2 = -(2/3) * (I2 ⊗ IC) + (2 * trC / 3) * ICodIC
    ℂ_dev = dγβ2 ⊗ V + γβ * dV2

    ℂ = ℂ_vol + ℂ_dev
    return _convect_tangent(ℂ, S, F)
end

p_wave_modulus(::Gent, props) = props[1] + 4 * props[2] / 3
