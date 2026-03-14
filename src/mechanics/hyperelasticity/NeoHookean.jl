"""
$(TYPEDEF)
"""
struct NeoHookean <: AbstractHyperelasticModel{2, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::NeoHookean, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    return [elastic_props.κ, elastic_props.μ]
end

"""
``\\psi = \\frac{1}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
        + \\frac{1}{2}\\mu\\left(\\bar{I}_1 - 3\\right)``
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::NeoHookean,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    # unpack properties
    κ, μ = props[1], props[2]

    # kinematics
    I       = one(∇u)
    F       = ∇u + I
    J       = det(F)
    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    I_1_bar = tr(J_m_23 * tdot(F))

    # constitutive
    ψ_vol = 0.5 * κ * (0.5 * (J * J - 1) - log(J))
    ψ_dev = 0.5 * μ * (I_1_bar - 3.)
    ψ     = ψ_vol + ψ_dev
    return ψ
end

function pk1_stress(
    ::NeoHookean,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    κ, μ    = props[1], props[2]
    F       = ∇u + one(∇u)
    J       = det(F)
    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    I_1     = tr(tdot(F))
    F_inv_T = inv(F)'
    P       = 0.5 * κ * (J * J - 1.) * F_inv_T +
              μ * J_m_23 * (F - (1. / 3.) * I_1 * F_inv_T)
    return P
end

# Analytical material tangent ∂P/∂F.
#
# Strategy: assemble the Lagrangian moduli ℂ = 2 ∂S/∂C (Norma convention) in
# terms of C⁻¹ outer products, then push forward via _convect_tangent.
#
# Volumetric:  ℂ_vol = κ [ J²(C⁻¹⊗C⁻¹) − (J²−1)(C⁻¹⊙C⁻¹) ]
# Deviatoric:  ℂ_dev = (2μ J⁻²/³/3) [ tr(C)(C⁻¹⊗C⁻¹/3 + C⁻¹⊙C⁻¹)
#                                       − I⊗C⁻¹ − C⁻¹⊗I ]
# where ⊙ is the minor-symmetrized outer product
#   (A⊙B)_{ijkl} = ½(A_{ik}B_{jl} + A_{il}B_{jk})
function material_tangent(
    ::NeoHookean,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    κ, μ  = props[1], props[2]
    F     = ∇u + one(∇u)
    J     = det(F)
    J2    = J * J
    C     = tdot(F)                          # F'F  (SymmetricTensor{2,3})
    IC    = inv(C)                           # C⁻¹  (SymmetricTensor{2,3})
    I2    = one(IC)                          # identity (SymmetricTensor{2,3})
    trC   = tr(C)
    Jm23  = cbrt(J)^(-2)

    # PK2 stress (needed for geometric stiffness in _convect_tangent)
    S = (0.5κ * (J2 - 1.0)) * IC + μ * Jm23 * (I2 - (trC / 3.0) * IC)

    # Lagrangian moduli  ℂ = 2 ∂S/∂C
    ICxIC  = otimes(IC, IC)                          # C⁻¹ ⊗ C⁻¹
    ICodIC = 0.5 * (otimesu(IC, IC) + otimesl(IC, IC))  # C⁻¹ ⊙ C⁻¹
    ℂ_vol  = κ * (J2 * ICxIC - (J2 - 1.0) * ICodIC)

    coeff  = 2.0μ * Jm23 / 3.0
    ℂ_dev  = coeff * (trC * (ICxIC / 3.0 + ICodIC) - otimes(IC, I2) - otimes(I2, IC))

    ℂ = ℂ_vol + ℂ_dev

    return _convect_tangent(ℂ, S, F)
end
