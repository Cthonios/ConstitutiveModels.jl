# Finite-deformation J2 (von Mises) plasticity with multiplicative split.
#
# Formulation
# -----------
#   Kinematics:    F = Fᵉ Fᵖ  (multiplicative decomposition)
#   Elasticity:    Hencky (logarithmic) stored energy in Mandel stress space
#   Yield surface: von Mises  f = σ_vm − (σ_y + H εᵖ)
#   Flow rule:     associated, isochoric  (tr Ṅ = 0)
#   Hardening:     linear isotropic
#
# State variables (NS = 10)
#   Z[1:9]  = vec(Fᵖ)  column-major; initial value = vec(I₃)
#   Z[10]   = εᵖ       accumulated equivalent plastic strain; initial 0
#
# Properties (NP = 4): [λ, μ, σ_y, H]
#   λ, μ  : Lamé constants for Hencky elastic response
#   σ_y   : initial yield stress
#   H     : linear isotropic hardening modulus
#
# GPU-compatible implementation
# ------------------------------
# log(Cᵉ_tr) and exp(Δεᵖ N) are computed via CM's branchless analytical
# 3×3 symmetric eigensolver (eigen_sym33_unit), applying scalar log/exp to
# the three eigenvalues.  N is symmetric and shares eigenvectors with Cᵉ_tr
# (both are isotropic functions of Cᵉ_tr), so the same spectral basis serves
# for both.  No LAPACK / LinearAlgebra calls remain.
#
# NOTE on state mutability
# ------------------------
# CM passes state as a mutable vector. Carina/FEC must allocate mutable
# per-quadrature-point state arrays of size NS=10 for this model.

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

# ---------------------------------------------------------------------------
# Internal: radial-return stress update
# ---------------------------------------------------------------------------
#
# Returns (W, P, state_new) where:
#   W         = stored energy density
#   P         = PK1 stress  (Tensor{2,3,T,9})
#   state_new = updated state (SVector{10,T})

@inline function _j2_stress(
    props,
    F::Tensor{2,3,T,9},
    state_old::AbstractVector,
) where T
    λ = T(props[1]); μ = T(props[2]); σ_y = T(props[3]); H = T(props[4])

    Fp_old   = Tensor{2,3,T,9}(ntuple(i -> T(state_old[i]), Val(9)))
    eqps     = T(state_old[10])

    # Trial elastic deformation gradient and right Cauchy-Green tensor
    Fe_tr    = F ⋅ inv(Fp_old)
    Ce_tr    = symmetric(tdot(Fe_tr))              # SymmetricTensor{2,3,T}

    # Hencky trial strain: Eᵉ_tr = ½ log(Cᵉ_tr)
    # eigen_sym33_unit is a GPU-compatible branchless analytical eigensolver
    Ee_tr    = _matrix_function(log, Ce_tr) / 2

    # Trial Mandel stress
    trEe     = tr(Ee_tr)
    I2       = one(SymmetricTensor{2,3,T})
    M_tr     = λ * trEe * I2 + 2μ * Ee_tr
    Mdev_tr  = dev(M_tr)
    σvm_tr   = sqrt(T(3) / 2) * norm(Mdev_tr)

    f_tr = σvm_tr - (σ_y + H * eqps)

    if f_tr ≤ zero(T)
        # Elastic step
        W        = T(0.5) * λ * trEe^2 + μ * dcontract(Ee_tr, Ee_tr)
        Fe_new   = Fe_tr
        Fp_new   = Fp_old
        eqps_new = eqps
        M_new    = M_tr
    else
        # Plastic step — radial return (exact for linear hardening)
        Δεᵖ      = f_tr / (3μ + H)
        N        = (T(3) / 2) * Mdev_tr / σvm_tr   # symmetric, traceless
        Ee_new   = Ee_tr - Δεᵖ * N
        eqps_new = eqps + Δεᵖ
        W        = T(0.5) * λ * trEe^2 + μ * dcontract(Ee_new, Ee_new)
        M_new    = λ * trEe * I2 + 2μ * Ee_new
        # exp(Δεᵖ N): N symmetric → spectral decomposition (GPU-compatible)
        exp_N    = Tensor{2,3,T,9}(_matrix_function(x -> exp(Δεᵖ * x), N))
        Fp_new   = exp_N ⋅ Fp_old
        Fe_new   = F ⋅ inv(Fp_new)
    end

    # PK1: P = Fᵉ⁻ᵀ M Fᵖ⁻ᵀ
    P = inv(Fe_new)' ⋅ Tensor{2,3,T,9}(M_new) ⋅ inv(Fp_new)'

    fp        = Fp_new.data
    state_new = SVector{10,T}(
        fp[1], fp[2], fp[3], fp[4], fp[5], fp[6], fp[7], fp[8], fp[9], eqps_new
    )
    return W, P, state_new
end

# ---------------------------------------------------------------------------
# Internal: analytical consistent tangent  ∂P/∂F
# ---------------------------------------------------------------------------
#
# Approximation: Fᵖ_new is FROZEN when differentiating P w.r.t. F.
#   • Elastic step (Fᵖ_new = Fᵖ_old, truly constant): result is exact.
#   • Plastic step: omitting ∂Fᵖ_new/∂F introduces O(Δεᵖ) relative error,
#     which is small in practice and preserves quadratic Newton convergence.
#
# ∂P/∂F = geometric + material contributions:
#
#   Geometric:  −(F⁻ᵀ)_{iL} P_{kJ}
#   Material:   2 Σ_{ABCD} (Fᵉ_new⁻ᵀ)_{iA} ℂ_{ABCD} (Fᵖ_new⁻ᵀ)_{BJ}
#                          (Fᵉ_tr)_{kC} (Fᵖ_old⁻¹)_{DL}
#
# ℂ assembled in the spectral basis {n_α} of Cᵉ_tr:
#
#   Material part:  Σ_{αβ} [ĉ_{αβ}/(2c_β)] (n_α⊗n_α)_{AB} (n_β⊗n_β)_{CD}
#   Geometric part: Σ_{α≠β} θ_{αβ} (n_α⊗n_β)_sym_{AB} (n_α⊗n_β)_sym_{CD}
#
# Principal-space algorithmic moduli ĉ_{αβ} = ∂m_α_new/∂ε_β_tr:
#   Elastic: ĉ_{αβ} = λ + 2μ δ_{αβ}
#   Plastic: ĉ_{αβ} = (λ + μβ̄) + 2μ(1 − 3β̄/2) δ_{αβ}
#                     + (2μβ̄ − 4μ²/(3μ+H)) N̄_α N̄_β
#            β̄ = 2μ Δεᵖ / σvm_tr,  N̄_α = 1.5 m_dev_α_tr / σvm_tr
#
# θ_{αβ} = (m_α_new − m_β_new)/(c_α − c_β)   [L'Hôpital if c_α ≈ c_β]

@inline function _j2_tangent_analytical(
    props,
    F::Tensor{2,3,T,9},
    state_old::AbstractVector,
    P::Tensor{2,3,T,9},
    state_new::AbstractVector,
) where T
    λ = T(props[1]); μ = T(props[2]); σ_y = T(props[3]); H = T(props[4])

    Fp_old    = Tensor{2,3,T,9}(ntuple(i -> T(state_old[i]), Val(9)))
    eqps      = T(state_old[10])
    Fp_new    = Tensor{2,3,T,9}(ntuple(i -> T(state_new[i]), Val(9)))
    Δεᵖ       = T(state_new[10]) - eqps

    Fp_old_i  = inv(Fp_old)
    Fp_new_i  = inv(Fp_new)
    Fe_tr     = F ⋅ Fp_old_i
    Fe_new    = F ⋅ Fp_new_i
    Fe_new_iT = inv(Fe_new)'
    F_iT      = inv(F)'

    # Spectral decomposition of symmetric trial Ce (GPU-compatible)
    Ce_tr = symmetric(tdot(Fe_tr))
    c, Q  = eigen_sym33_unit(Ce_tr)   # c::Vec{3,T}, Q::Tensor{2,3,T,9}

    # Trial principal Hencky strains and Mandel stresses
    ε_tr    = SVector{3,T}(T(0.5)*log(c[1]), T(0.5)*log(c[2]), T(0.5)*log(c[3]))
    tr_ε    = ε_tr[1] + ε_tr[2] + ε_tr[3]
    m_tr    = SVector{3,T}(
        λ*tr_ε + 2μ*ε_tr[1],
        λ*tr_ε + 2μ*ε_tr[2],
        λ*tr_ε + 2μ*ε_tr[3],
    )
    tr_m    = m_tr[1] + m_tr[2] + m_tr[3]
    mdev_tr = SVector{3,T}(m_tr[1] - tr_m/3, m_tr[2] - tr_m/3, m_tr[3] - tr_m/3)
    σvm_tr  = sqrt(T(3)/2 * (mdev_tr[1]^2 + mdev_tr[2]^2 + mdev_tr[3]^2))

    f_tr = σvm_tr - (σ_y + H * eqps)

    if f_tr ≤ zero(T)
        # Elastic algorithmic moduli: ĉ[α,β] = λ + 2μ δ(α,β)
        m_new = m_tr
        ĉ = SMatrix{3,3,T,9}(λ+2μ, λ, λ,  λ, λ+2μ, λ,  λ, λ, λ+2μ)
    else
        # Plastic algorithmic moduli
        N̄ = SVector{3,T}(
            T(1.5)*mdev_tr[1]/σvm_tr,
            T(1.5)*mdev_tr[2]/σvm_tr,
            T(1.5)*mdev_tr[3]/σvm_tr,
        )
        m_new = SVector{3,T}(
            m_tr[1] - 2μ*Δεᵖ*N̄[1],
            m_tr[2] - 2μ*Δεᵖ*N̄[2],
            m_tr[3] - 2μ*Δεᵖ*N̄[3],
        )
        β̄     = 2μ * Δεᵖ / σvm_tr
        Ac    = λ + μ * β̄
        Bc    = 2μ * (1 - T(1.5) * β̄)
        Ccoef = 2μ * β̄ - 4μ^2 / (3μ + H)
        ĉ = SMatrix{3,3,T,9}(ntuple(Val(9)) do lin
            i = (lin - 1) % 3 + 1
            j = (lin - 1) ÷ 3 + 1
            Ac + Bc * ifelse(i == j, one(T), zero(T)) + Ccoef * N̄[i] * N̄[j]
        end)
    end

    # Assemble ℂ_{ABCD} in principal basis
    CC = MArray{Tuple{3,3,3,3},T,4,81}(ntuple(_ -> zero(T), Val(81)))

    # Material part: Σ_{αβ} [ĉ_{αβ}/(2c_β)] (nα⊗nα)_{AB} (nβ⊗nβ)_{CD}
    for α in 1:3, β in 1:3
        coeff = ĉ[α, β] / (2 * c[β])
        nα = SVector{3,T}(Q[1,α], Q[2,α], Q[3,α])
        nβ = SVector{3,T}(Q[1,β], Q[2,β], Q[3,β])
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            CC[A, B, C, D] += coeff * nα[A] * nα[B] * nβ[C] * nβ[D]
        end
    end

    # Geometric part (α ≠ β): θ_{αβ} (nα⊗nβ)sym_{AB} (nα⊗nβ)sym_{CD}
    tol = T(1e-10)
    for α in 1:3, β in 1:3
        α == β && continue
        nα = SVector{3,T}(Q[1,α], Q[2,α], Q[3,α])
        nβ = SVector{3,T}(Q[1,β], Q[2,β], Q[3,β])
        Δc = c[α] - c[β]
        θ = if abs(Δc) > tol * (c[α] + c[β])
            (m_new[α] - m_new[β]) / Δc
        else
            # L'Hôpital limit
            (ĉ[α, α] - ĉ[β, α]) / (2 * c[α])
        end
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            sAB = T(0.5) * (nα[A] * nβ[B] + nα[B] * nβ[A])
            sCD = T(0.5) * (nα[C] * nβ[D] + nα[D] * nβ[C])
            CC[A, B, C, D] += θ * sAB * sCD
        end
    end

    # Assemble AA_{iJkL} = −(F⁻ᵀ)_{iL} P_{kJ}
    #   + 2 Σ_{ABCD} (Fᵉ_new⁻ᵀ)_{iA} ℂ_{ABCD} (Fᵖ_new⁻ᵀ)_{BJ} (Fᵉ_tr)_{kC} (Fᵖ_old⁻¹)_{DL}
    Fp_new_iT = Fp_new_i'
    AA = MArray{Tuple{3,3,3,3},T,4,81}(undef)
    for i in 1:3, J in 1:3, k in 1:3, L in 1:3
        val = -F_iT[i, L] * P[k, J]
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            val += 2 * Fe_new_iT[i, A] * CC[A, B, C, D] *
                   Fp_new_iT[B, J] * Fe_tr[k, C] * Fp_old_i[D, L]
        end
        AA[i, J, k, L] = val
    end
    return Tensor{4,3,T,81}(ntuple(i -> AA[i], Val(81)))
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
    W, _, state_new_vec = _j2_stress(props, F, Z_old)
    Z_new .= state_new_vec
    return W
end

function pk1_stress(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    _, P, state_new_vec = _j2_stress(props, F, Z_old)
    Z_new .= state_new_vec
    return P
end

function material_tangent(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F = ∇u + one(∇u)
    _, P, state_new_vec = _j2_stress(props, F, Z_old)
    Z_new .= state_new_vec
    return _j2_tangent_analytical(props, F, Z_old, P, state_new_vec)
end
