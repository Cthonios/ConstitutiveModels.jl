struct ArrudaBoyce <: AbstractHyperelasticModel{3, 0}
end

function initialize_props(::ArrudaBoyce, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    n = inputs["n"]
    return [elastic_props.κ, elastic_props.μ, sqrt(n)]
end

# using treloar approximation
# function _inverse_langevin_approx(y)
#     return 3 * y / (1 - (3 / 5 * y^2 + 36 / 175 * y^4 + 108 / 875 * y^6))
# end

# function _dinverse_langevin_approx(y)
#     return 13125 * (108 * y^6 + 108 * y^4 + 105 * y^2 + 175) /
#         (-108 * y^6 - 180 * y^4 - 525 * y^2 + 875)^2
# end

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

# function pk1_stress(
#     ::ArrudaBoyce,
#     props, Δt,
#     ∇u, θ, Z_old, Z_new
# )
#     κ, μ, sqrt_n = props[1], props[2], props[3]

#     # kinematics
#     I       = one(typeof(∇u))
#     F       = ∇u + I
#     J       = det(F)
#     J_m_13  = 1. / cbrt(J)
#     J_m_23  = J_m_13 * J_m_13
#     C       = tdot(F)
#     I1      = tr(C)
#     F_inv_T = inv(F)'

#     # chain stretch
#     I1_bar  = J_m_23 * I1
#     λ_chain = sqrt(I1_bar / 3)

#     β = inverse_langevin_approximation(
#         λ_chain / sqrt_n,
#         TreloarApproximation()
#     )

#     # volumetric part
#     P_vol = 0.5 * κ * (J * J - 1.) * F_inv_T

#     # deviatoric part
#     coeff = 2μ * β / (3 * λ_chain)

#     P_dev = coeff * J_m_23 * (
#         F -
#         (1 / 3) * I1 * F_inv_T
#     )

#     P = P_vol + P_dev
#     return P
# end
