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
# Properties (NP = 5): [ρ, λ, μ, σ_y, H]
#   ρ     : Lagrangian-frame density; every model carries it as props[1]
#   λ, μ  : Lamé constants (converted to κ = λ + 2μ/3 internally)
#   σ_y   : initial yield stress
#   H     : linear isotropic hardening modulus

"""
$(TYPEDEF)
"""
struct FiniteDefJ2Plasticity <: AbstractConstitutiveModel
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::FiniteDefJ2Plasticity, inputs::Dict{String})
    ρ   = get_property(inputs, "density")
    ec  = ElasticConstants(inputs)
    σ_y = get_property(inputs, "yield stress", 0.0)
    H   = get_property(inputs, "hardening modulus", 0.0)
    return [ρ, ec.λ, ec.μ, σ_y, H]
end

"""
Initial state: Fᵖ = I₃, εᵖ = 0.
$(TYPEDSIGNATURES)
"""
function initialize_state(::FiniteDefJ2Plasticity, float_type = Float64)
    return float_type[1, 0, 0, 0, 1, 0, 0, 0, 1, 0]   # vec(I) ++ 0
end

num_properties(::FiniteDefJ2Plasticity) = 5
num_state_variables(::FiniteDefJ2Plasticity) = 10

function property_names(::FiniteDefJ2Plasticity)
    return [
        "density", "Lamé's first constant", "shear modulus",
        "yield stress", "hardening modulus"
    ]
end

function state_variable_names(::FiniteDefJ2Plasticity)
    return [
        "Fp_xx", "Fp_yx", "Fp_zx",
        "Fp_xy", "Fp_yy", "Fp_zy",
        "Fp_xz", "Fp_yz", "Fp_zz",
        "eqps",
    ]
end

# ---------------------------------------------------------------------------
# Internal: Simo-Hughes stress update (BOX 9.1)
# ---------------------------------------------------------------------------

# `κ` is the bulk modulus used for the volumetric response; the public
# functions pass the model's own κ = λ + 2μ/3, and the `isochoric_*` functions
# pass zero, which removes the volumetric energy and stress and leaves the
# return map unchanged (it acts on the isochoric b̄ᵉ alone).
@inline function _sh_j2_stress(
    props,
    F::Tensor{2,3,T,9},
    state_old::AbstractVector,
    κ::T,
) where T
    μ = T(props[3]); σ_y = T(props[4]); K = T(props[5])

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
#
# What this returns is the major-symmetric part of the Jacobian of the stress
# update, not the Jacobian itself.  The return map evaluates the effective
# shear modulus μ̄ = μ tr(b̄ᵉ_trial)/3 at the trial state, so the discrete
# update is not exactly variational and its Jacobian carries an antisymmetric
# component.  BOX 9.2 drops that component, which is what a solver assembling
# a symmetric operator wants.
#
# The coefficient β₂ of BOX 9.2 is (1 - 1/β₀) (2/3) ‖s_trial‖ Δγ / μ̄, with
# the effective modulus μ̄, not μ.  An earlier version divided by μ; the two
# agree when tr(b̄ᵉ_trial) = 3, so the error was invisible at small strain and
# with H = 0 (β₀ = 1), and grew with the hardening modulus and with
# tr(b̄ᵉ_trial) - 3: 1.7e-5 relative in the symmetric part at an equivalent
# plastic strain of 0.11 (one increment from the virgin state, E = 1e3,
# ν = 0.3, σ_y = 20, H = 100) and 9.8e-2 at 0.63 (tr(b̄ᵉ_trial) - 3 = 3.6).
#
# Measured against a central-difference Jacobian of `pk1_stress` (h = 1e-6),
# as ‖·‖ relative to the tangent norm, at the same material and increments:
#
#        H     eqps    ‖sym(A_fd) - A‖/‖A‖    ‖asym(A_fd)‖/‖A_fd‖
#       100    0.11          3.6e-11               3.7e-03
#       100    0.37          1.0e-10               4.7e-03
#       100    0.63          1.4e-10               5.3e-03
#       300    0.59          1.1e-10               3.8e-03
#
# The symmetric part is recovered to the accuracy of the differences and the
# remaining discrepancy is the antisymmetric part; the elastic branch is the
# exact Jacobian.  Newton therefore stays fast on steps that yield but is not
# fully quadratic there.

# One step of the pull-back: contract index `k` of the fourth-order tensor X
# with F⁻¹, X'[..., K, ...] = Σ_x F⁻¹[K, x] X[..., x, ...].  Built from an
# ntuple into an immutable Tensor, as _convect_tangent is: a mutable 81-entry
# buffer (MArray) is allocated on the device heap by GPU compilers, once per
# call, which exhausts the CUDA malloc heap in a tangent kernel.
@inline function _pull_back_index(F_inv::Tensor{2,3,T,9}, X, ::Val{k}) where {T, k}
    data = ntuple(Val(81)) do lin
        l, rem = divrem(lin - 1, 27)
        kk, rem = divrem(rem, 9)
        j, i = divrem(rem, 3)
        idx = (i + 1, j + 1, kk + 1, l + 1)
        v = zero(T)
        for x in 1:3
            jdx = Base.setindex(idx, x, k)
            v += F_inv[idx[k], x] * X[jdx[1], jdx[2], jdx[3], jdx[4]]
        end
        v
    end
    return Tensor{4, 3, T, 81}(data)
end

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
    κ::T,
) where T
    μ = T(props[3]); σ_y = T(props[4]); K = T(props[5])

    J = det(F)
    F_inv = inv(F)
    S = F_inv ⋅ P   # 2nd Piola-Kirchhoff

    I2 = one(SymmetricTensor{2, 3, T, 6})
    I4_sym = one(SymmetricTensor{4, 3, T, 36})

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
        β₂ = (one(T) - one(T)/β₀) * T(2)/3 * s_trial_norm / μ̄ * Δγ
        β₃ = one(T)/β₀ - β₁ + β₂
        β₄ = (one(T)/β₀ - β₁) * s_trial_norm / μ̄

        n_sq = symmetric(n ⋅ n)
        c_dev_n2 = symmetric(n ⊗ dev(n_sq) + dev(n_sq) ⊗ n) / 2

        CC_spatial = c_vol + c_dev_trial -
                     β₁ * c_dev_trial -
                     2μ̄ * β₃ * (n ⊗ n) -
                     2μ̄ * β₄ * c_dev_n2
    end

    # Pull-back: spatial → material,
    #   CC[A,B,C,D] = F⁻¹[A,a] F⁻¹[B,b] c[a,b,c,d] F⁻¹[C,c] F⁻¹[D,d]
    #
    # Contracting one index at a time costs 4·3⁵ multiplies; contracting all
    # four at once costs 3⁸ -- about thirty times more for the same result.
    T1 = _pull_back_index(F_inv, CC_spatial, Val(1))
    T2 = _pull_back_index(F_inv, T1, Val(2))
    T3 = _pull_back_index(F_inv, T2, Val(3))
    CC = _pull_back_index(F_inv, T3, Val(4))

    return _convect_tangent(CC, S, F)
end

# ---------------------------------------------------------------------------
# CM public API
# ---------------------------------------------------------------------------

# Bulk modulus κ = λ + 2μ/3 from the property vector [ρ, λ, μ, σ_y, H].
@inline _j2_bulk_modulus(props, ::Type{T}) where T = T(props[2]) + 2 * T(props[3]) / 3

# The state container is a view into the assembler's storage with a length
# known only at run time.  Broadcasting an SVector into it (`Z_new .= vec`)
# carries the DimensionMismatch path of the broadcast shape check, whose
# message is built with `show`; that path cannot be compiled for a GPU, so
# the entries are copied one by one.
@inline function _store_state!(Z_new, state_new_vec::SVector{N, T}) where {N, T}
    for i in 1:N
        Z_new[i] = state_new_vec[i]
    end
    return nothing
end

@inline function _j2_energy(props, Z_old, Z_new, ∇u, κ)
    F = ∇u + one(∇u)
    W, _, state_new_vec, _, _, _, _, _, _ = _sh_j2_stress(props, F, Z_old, κ)
    _store_state!(Z_new, state_new_vec)
    return W
end

@inline function _j2_pk1(props, Z_old, Z_new, ∇u, κ)
    F = ∇u + one(∇u)
    _, P, state_new_vec, _, _, _, _, _, _ = _sh_j2_stress(props, F, Z_old, κ)
    _store_state!(Z_new, state_new_vec)
    return P
end

# The stress and the tangent from one return map.
@inline function _j2_stress_tangent(props, Z_old, Z_new, ∇u, κ)
    F = ∇u + one(∇u)
    W, P, state_new_vec, s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n =
        _sh_j2_stress(props, F, Z_old, κ)
    _store_state!(Z_new, state_new_vec)
    A = _sh_j2_tangent(props, F, Z_old, P,
                       s_new, be_bar_tr, s_trial_norm, μ̄, Δγ, α_n, κ)
    return P, A
end

@inline _j2_tangent(props, Z_old, Z_new, ∇u, κ) = _j2_stress_tangent(props, Z_old, Z_new, ∇u, κ)[2]

function helmholtz_free_energy(
    ::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    return _j2_energy(props, Z_old, Z_new, ∇u, _j2_bulk_modulus(props, eltype(∇u)))
end

function pk1_stress(
    ::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    return _j2_pk1(props, Z_old, Z_new, ∇u, _j2_bulk_modulus(props, eltype(∇u)))
end

function material_tangent(
    ::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    return _j2_tangent(props, Z_old, Z_new, ∇u, _j2_bulk_modulus(props, eltype(∇u)))
end

function pk1_stress_and_material_tangent(
    ::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    return _j2_stress_tangent(props, Z_old, Z_new, ∇u, _j2_bulk_modulus(props, eltype(∇u)))
end

# ---------------------------------------------------------------------------
# Volumetric-isochoric split (see Interface.jl)
#
#   W = κ/2 (J - 1)² + μ/2 (tr b̄ᵉ - 3),   θ(J) = J - 1,   θ'(J) = 1,
#
# the split of BOX 9.1: det b̄ᵉ = 1 makes the isochoric term independent of J.
# ---------------------------------------------------------------------------

has_volumetric_isochoric_split(::FiniteDefJ2Plasticity) = true
volumetric_strain(::FiniteDefJ2Plasticity, J) = J - one(J)
volumetric_strain_derivative(::FiniteDefJ2Plasticity, J) = one(J)
volumetric_strain_second_derivative(::FiniteDefJ2Plasticity, J) = zero(J)
bulk_modulus(::FiniteDefJ2Plasticity, props) = _j2_bulk_modulus(props, eltype(props))

function isochoric_helmholtz_free_energy(
    ::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    return _j2_energy(props, Z_old, Z_new, ∇u, zero(eltype(∇u)))
end

function isochoric_pk1_stress(
    ::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    return _j2_pk1(props, Z_old, Z_new, ∇u, zero(eltype(∇u)))
end

function isochoric_material_tangent(
    ::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    return _j2_tangent(props, Z_old, Z_new, ∇u, zero(eltype(∇u)))
end

"""
Cauchy stress σ = J⁻¹ P Fᵀ.

`FiniteDefJ2Plasticity` subtypes `AbstractConstitutiveModel` directly rather
than `AbstractHyperelasticModel` -- it is path dependent, so the hyperelastic
defaults (AD tangents from a stored energy) do not apply to it.  That also means
it does not inherit the hyperelastic `cauchy_stress` fallback, so the push-forward
is written out here.  Consumers that report stress (e.g. Carina's output writer)
call this.
$(TYPEDSIGNATURES)
"""
function cauchy_stress(
    model::FiniteDefJ2Plasticity,
    props, Z_old, Z_new, Δt,
    ∇u, θ
)
    F = ∇u + one(∇u)
    J = det(F)
    P = pk1_stress(model, props, Z_old, Z_new, Δt, ∇u, θ)
    return (1 / J) * dot(P, transpose(F))
end

p_wave_modulus(::FiniteDefJ2Plasticity, props) = props[2] + 2 * props[3]
