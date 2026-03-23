# Simo-Hughes J2 finite-deformation plasticity (BOX 9.1 + 9.2).
#
# Reference: Simo & Hughes, Computational Inelasticity, pp 317-321.
#
# Formulation
# -----------
#   Kinematics:    F = Fᵉ Fᵖ  (multiplicative decomposition)
#   Elasticity:    neo-Hookean-type vol/dev split: U(J) + μ/2(tr b̄ᵉ - 3)
#   Yield surface: von Mises  f = ‖s‖ - √(2/3)(σ_y + K α)
#   Flow rule:     associated, isochoric  (tr Ṅ = 0)
#   Hardening:     linear isotropic
#
# State variables (NS = 10)
#   Z[1:9]  = vec(Fᵖ)  column-major; initial value = vec(I₃)
#   Z[10]   = α        accumulated equivalent plastic strain; initial 0
#
# Properties (NP = 4): [λ, μ, σ_y, H]
#   λ, μ  : Lamé constants (converted to κ = λ + 2μ/3 internally)
#   σ_y   : initial yield stress
#   H     : linear isotropic hardening modulus

"""
$(TYPEDEF)
"""
struct FiniteDefJ2Plasticity <: AbstractMechanicalModel{4, 10}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::FiniteDefJ2Plasticity, inputs::Dict{String})
    ec  = ElasticConstants(inputs)
    σ_y = get(inputs, "yield stress",       0.0)
    H   = get(inputs, "hardening modulus",  0.0)
    return [ec.λ, ec.μ, σ_y, H]
end

"""
Initial state: Fᵖ = I₃, εᵖ = 0.
$(TYPEDSIGNATURES)
"""
function initialize_state(::FiniteDefJ2Plasticity, float_type = Float64)
    return float_type[1, 0, 0, 0, 1, 0, 0, 0, 1, 0]   # vec(I) ++ 0
end

function state_variable_names(::FiniteDefJ2Plasticity)
    return [
        "Fp_11", "Fp_21", "Fp_31",
        "Fp_12", "Fp_22", "Fp_32",
        "Fp_13", "Fp_23", "Fp_33",
        "eqps",
    ]
end

# ---------------------------------------------------------------------------
# Internal: Simo-Hughes stress update (BOX 9.1)
# ---------------------------------------------------------------------------

@inline function _sh_j2_stress(
    props,
    F::Tensor{2,3,T,9},
    state_old::AbstractVector,
) where T
    # Convert Lamé λ → bulk modulus κ = λ + 2μ/3
    λ = T(props[1]); μ = T(props[2]); σ_y = T(props[3]); K = T(props[4])
    κ = λ + 2μ / 3

    Fp_old = Tensor{2,3,T,9}(ntuple(i -> T(state_old[i]), Val(9)))
    α_n    = T(state_old[10])

    # Total Jacobian and isochoric factor
    J    = det(F)
    Jm23 = J^(-T(2)/3)

    # Trial elastic left Cauchy-Green (isochoric): b̄ᵉ_trial = J^{-2/3} Fe_tr · Fe_trᵀ
    Fe_tr     = F ⋅ inv(Fp_old)
    be_bar_tr = symmetric(Jm23 * (Fe_tr ⋅ Fe_tr'))

    # Trial deviatoric Kirchhoff stress: s_trial = μ dev[b̄ᵉ_trial]
    s_trial      = μ * dev(be_bar_tr)
    s_trial_norm = norm(s_trial)

    # Effective shear modulus: μ̄ = μ/3 tr[b̄ᵉ_trial]
    μ̄ = μ * tr(be_bar_tr) / 3

    # Yield function
    f_trial = s_trial_norm - sqrt(T(2)/3) * (σ_y + K * α_n)

    I2 = one(SymmetricTensor{2,3,T})

    if f_trial ≤ zero(T)
        # Elastic step
        s_new      = s_trial
        be_bar_new = be_bar_tr
        α_new      = α_n
        Δγ         = zero(T)
        Fp_new     = Fp_old
    else
        # Plastic step — radial return (BOX 9.1, step 4)
        n  = s_trial / s_trial_norm
        Δγ = f_trial / (2μ̄ + T(2)/3 * K)

        s_new = s_trial - 2μ̄ * Δγ * n
        α_new = α_n + sqrt(T(2)/3) * Δγ

        # Update b̄ᵉ (eq 9.3.33)
        Ie_bar     = tr(be_bar_tr) / 3
        be_bar_new = s_new / μ + Ie_bar * I2

        # Recover Fp_new from b̄ᵉ_new via polar decomposition
        be_tr_sqrt_inv = Tensor{2,3,T,9}(_matrix_function(x -> 1/sqrt(x), be_bar_tr))
        Fe_tr_iso      = J^(-T(1)/3) * Fe_tr
        R_tr           = be_tr_sqrt_inv ⋅ Fe_tr_iso

        be_new_sqrt    = Tensor{2,3,T,9}(_matrix_function(sqrt, be_bar_new))
        Fe_new_iso     = be_new_sqrt ⋅ R_tr
        Fe_new         = J^(T(1)/3) * Fe_new_iso
        Fp_new         = inv(Fe_new) ⋅ F
    end

    # Kirchhoff stress: τ = J p I + s,  p = κ(J-1)
    p = κ * (J - one(T))
    τ = J * p * I2 + s_new

    # PK1: P = τ · F⁻ᵀ
    P = Tensor{2,3,T,9}(τ) ⋅ inv(F)'

    # Energy
    W = κ / 2 * (J - 1)^2 + μ / 2 * (tr(be_bar_new) - T(3))

    fp = Fp_new.data
    state_new = SVector{10,T}(
        fp[1], fp[2], fp[3], fp[4], fp[5], fp[6], fp[7], fp[8], fp[9], α_new
    )
    return W, P, state_new, s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n
end

# ---------------------------------------------------------------------------
# Internal: Simo-Hughes consistent tangent (BOX 9.2)
# Uses _convect_tangent from utils/TensorUtils.jl for push-forward.
# ---------------------------------------------------------------------------

@inline function _sh_j2_tangent(
    props,
    F::Tensor{2,3,T,9},
    state_old::AbstractVector,
    P::Tensor{2,3,T,9},
    s_new::SymmetricTensor{2,3,T},
    be_bar_tr::SymmetricTensor{2,3,T},
    s_trial_norm::T,
    μ̄::T,
    Δγ::T,
    α_n::T,
) where T
    λ = T(props[1]); μ = T(props[2]); σ_y = T(props[3]); K = T(props[4])
    κ = λ + 2μ / 3

    J = det(F)
    F_inv = inv(F)
    S = F_inv ⋅ P   # 2nd Piola-Kirchhoff

    I2 = one(SymmetricTensor{2,3,T})
    I4_sym = one(SymmetricTensor{4,3,T})

    # Volumetric spatial tangent
    coeff_1x1 = κ * J * (2*J - one(T))
    coeff_I   = 2 * κ * J * (J - one(T))
    c_vol = coeff_1x1 * (I2 ⊗ I2) - coeff_I * I4_sym

    # Unit normal
    s_trial_recomp = μ * dev(be_bar_tr)
    n = s_trial_norm > zero(T) ? s_trial_recomp / s_trial_norm :
        zero(SymmetricTensor{2,3,T})

    # Deviatoric trial tangent
    c_dev_trial = 2μ̄ * (I4_sym - T(1)/3 * (I2 ⊗ I2)) -
                  T(2)/3 * s_trial_norm * (n ⊗ I2 + I2 ⊗ n)

    f_trial = s_trial_norm - sqrt(T(2)/3) * (σ_y + K * α_n)

    if f_trial ≤ zero(T)
        CC_spatial = c_vol + c_dev_trial
    else
        # Plastic correction (BOX 9.2, steps 2-3)
        β₀ = one(T) + K / (3μ̄)
        β₁ = 2μ̄ * Δγ / s_trial_norm
        β₂ = (one(T) - one(T)/β₀) * T(2)/3 * s_trial_norm / μ * Δγ
        β₃ = one(T)/β₀ - β₁ + β₂
        β₄ = (one(T)/β₀ - β₁) * s_trial_norm / μ̄

        n_sq = symmetric(n ⋅ n)
        c_dev_n2 = symmetric(n ⊗ dev(n_sq) + dev(n_sq) ⊗ n) / 2

        CC_spatial = c_vol + c_dev_trial -
                     β₁ * c_dev_trial -
                     2μ̄ * β₃ * (n ⊗ n) -
                     2μ̄ * β₄ * c_dev_n2
    end

    # Pull-back: spatial → material
    CC = MArray{Tuple{3,3,3,3},T,4,81}(ntuple(_ -> zero(T), Val(81)))
    for A in 1:3, B in 1:3, C in 1:3, D in 1:3
        val = zero(T)
        for a in 1:3, b in 1:3, c in 1:3, d in 1:3
            val += F_inv[A, a] * F_inv[B, b] * CC_spatial[a, b, c, d] * F_inv[C, c] * F_inv[D, d]
        end
        CC[A, B, C, D] = val
    end

    return _convect_tangent(CC, S, F)
end

# ---------------------------------------------------------------------------
# CM public API
# ---------------------------------------------------------------------------

function helmholtz_free_energy(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    W, _, state_new_vec, _, _, _, _, _, _ = _sh_j2_stress(props, F, Z_old)
    Z_new .= state_new_vec
    return W
end

function pk1_stress(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    _, P, state_new_vec, _, _, _, _, _, _ = _sh_j2_stress(props, F, Z_old)
    Z_new .= state_new_vec
    return P
end

function material_tangent(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    W, P, state_new_vec, s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n =
        _sh_j2_stress(props, F, Z_old)
    Z_new .= state_new_vec
    return _sh_j2_tangent(props, F, Z_old, P,
                           s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n)
end
