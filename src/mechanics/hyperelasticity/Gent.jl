"""
$(TYPEDEF)
"""
struct Gent <: AbstractHyperelasticModel{3, 0}
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
function helmholtz_free_energy(
    ::Gent,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    # unpack properties
    κ, μ, Jm = props[1], props[2], props[3]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
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

function pk1_stress(
    ::Gent,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    # unpack properties
    κ, μ, Jm = props[1], props[2], props[3]

    # kinematics
    F       = ∇u + one(typeof(∇u))
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
