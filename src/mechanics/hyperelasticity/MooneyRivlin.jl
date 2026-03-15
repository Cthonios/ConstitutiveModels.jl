"""
$(TYPEDEF)
"""
struct MooneyRivlin <: AbstractHyperelasticModel{3,0}  # 3 material props: κ, C1, C2
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::MooneyRivlin, inputs::Dict{String})
    ec = ElasticConstants(inputs)
    C1 = get(inputs, "C1", ec.μ / 2)  # default split of μ
    C2 = get(inputs, "C2", ec.μ / 2)
    return [ec.κ, C1, C2]
end

"""
``\\psi = \\frac{1}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
        + C_1\\left(\\bar{I}_1 - 3\\right)
        + C_2\\left(\\bar{I}_2 - 3\\right)
``
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::MooneyRivlin,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    κ, C1, C2 = props[1], props[2], props[3]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
    J       = det(F)
    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    C       = tdot(F)

    # invariants
    I1 = tr(C)
    I2 = 0.5 * (I1 * I1 - tr(C ⋅ C))
    I1_bar = J_m_23 * I1
    I2_bar = J_m_23 * J_m_23 * I2

    # energies
    ψ_vol = 0.5 * κ * (0.5 * (J^2 - 1) - log(J))
    ψ_dev = C1 * (I1_bar - 3.0) + C2 * (I2_bar - 3.0)
    return ψ_vol + ψ_dev
end

"""
$(TYPEDSIGNATURES)
"""
function pk1_stress(
    ::MooneyRivlin,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    κ, C1, C2 = props[1], props[2], props[3]

    # kinematics
    F       = ∇u + one(∇u)
    J       = det(F)
    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    J_m_43  = J_m_23 * J_m_23
    C       = tdot(F)
    I1      = tr(C)
    I2      = 0.5 * (I1^2 - tr(C ⋅ C))
    F_inv_T = inv(F)'

    # constitutive
    P_vol = 0.5 * κ * (J * J - 1.) * F_inv_T
    P_dev = C1 * J_m_23 * (
        2.0 * F - 
        (2.0 / 3.0) * I1 * F_inv_T
    ) + C2 * J_m_43 * (
        2.0 * (I1 * F - F ⋅ C) - 
        (4.0 / 3.0) * I2 * F_inv_T
    )
    P = P_vol + P_dev
    return P
end


# TODO fix this, failing on the I2 term
# things are passing when we set C2 = 0 since
# that's simply NeoHookean
# """
# $(TYPEDSIGNATURES)
# """
# function material_tangent(
#     ::MooneyRivlin,
#     props, Δt,
#     ∇u, θ, Z_old, Z_new
# )
#     κ, C1, C2 = props

#     # Kinematics
#     F    = ∇u + one(∇u)
#     J    = det(F)
#     J2   = J * J
#     C    = tdot(F)
#     IC   = inv(C)
#     I    = one(IC)
#     trC  = tr(C)
#     I2   = 0.5 * (trC^2 - tr(C ⋅ C))
#     Jm13 = 1 / cbrt(J)
#     Jm23 = Jm13^2
#     Jm43 = Jm23^2

#     # Pk2 stress
#     S = (0.5κ * (J2 - 1.0)) * IC + 
#         2C1 * Jm23 * (one(IC) - (trC / 3.0) * IC) + 
#         2C2 * Jm43 * (trC * I - C - (2.0 / 3.0) * I2 * IC)

#     # Lagrangian moduli  ℂ = 2 ∂S/∂C
#     ICxIC  = otimes(IC, IC)                          # C⁻¹ ⊗ C⁻¹
#     ICodIC = 0.5 * (otimesu(IC, IC) + otimesl(IC, IC))  # C⁻¹ ⊙ C⁻¹
#     ℂ_vol  = κ * (J2 * ICxIC - (J2 - 1.0) * ICodIC)

#     coeff1  = 4C1 * Jm23 / 3.0
#     coeff2 = 2C2 * Jm43
#     ℂ_dev  = coeff1 * (
#         trC * (ICxIC / 3.0 + ICodIC) - 
#         IC ⊗ I - 
#         I ⊗ IC
#     ) + coeff2 * (
#         trC * (ICxIC / 3.0 + ICodIC) -
#         ICodIC -
#         IC ⊗ C -
#         C ⊗ IC
#     )

#     ℂ = ℂ_vol + ℂ_dev

#     return _convect_tangent(ℂ, S, F)
# end
