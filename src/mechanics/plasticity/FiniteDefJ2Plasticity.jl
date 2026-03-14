import LinearAlgebra: I, Symmetric, eigen, norm, tr, exp, inv, log

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
# NOTE on GPU / static-array support
# ------------------------------------
# The stress update uses log(Cᵉ_tr) and exp(Δεᵖ N) for which GPU-compatible
# static-array implementations are not yet available. This model runs on CPU
# only. GPU support is a future extension requiring a static-array matrix
# exponential (e.g. Cayley-Hamilton closed form for traceless 3×3).
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
#   P         = PK1 stress  (SMatrix{3,3,Float64,9})
#   state_new = updated state vector (length 10)

function _j2_stress(
    props,
    F_sa::SMatrix{3,3,Float64,9},
    state_old::AbstractVector{Float64},
)
    λ, μ, σ_y, H = props[1], props[2], props[3], props[4]

    Fp_old  = Matrix{Float64}(reshape(state_old[1:9], 3, 3))
    eqps    = state_old[10]

    # Trial elastic deformation gradient
    Fe_tr   = Matrix{Float64}(F_sa) * inv(Fp_old)
    Ce_tr   = Symmetric(Fe_tr' * Fe_tr)
    Ee_tr   = 0.5 * log(Matrix(Ce_tr))      # Hencky trial strain

    # Trial Mandel stress
    trEe    = tr(Ee_tr)
    M_tr    = λ * trEe * I(3) + 2μ * Ee_tr
    Mdev_tr = M_tr - (tr(M_tr) / 3) * I(3)
    σvm_tr  = sqrt(1.5) * norm(Mdev_tr)

    f_tr = σvm_tr - (σ_y + H * eqps)

    if f_tr ≤ 0.0
        # Elastic step
        W       = 0.5λ * trEe^2 + μ * tr(Ee_tr * Ee_tr)
        Fe_new  = Fe_tr
        Fp_new  = Fp_old
        eqps_new = eqps
        M_new   = M_tr
    else
        # Plastic step — radial return (exact for linear hardening)
        Δεᵖ     = f_tr / (3μ + H)
        N       = 1.5 * Mdev_tr / σvm_tr     # tr(N) = 0
        Fp_new  = exp(Δεᵖ * N) * Fp_old
        eqps_new = eqps + Δεᵖ
        Ee_new  = Ee_tr - Δεᵖ * N
        W       = 0.5λ * trEe^2 + μ * tr(Ee_new * Ee_new)
        M_new   = λ * trEe * I(3) + 2μ * Ee_new
        Fe_new  = Matrix{Float64}(F_sa) * inv(Fp_new)
    end

    # PK1:  P = Fᵉ⁻ᵀ M Fᵖ⁻ᵀ
    P_dense = inv(Fe_new)' * M_new * inv(Fp_new)'
    P = SMatrix{3,3,Float64,9}(P_dense)

    state_new = Vector{Float64}(undef, 10)
    state_new[1:9] = vec(Fp_new)
    state_new[10]  = eqps_new

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

function _j2_tangent_analytical(
    props,
    F_sa::SMatrix{3,3,Float64,9},
    state_old::AbstractVector{Float64},
    P::SMatrix{3,3,Float64,9},
    state_new::AbstractVector{Float64},
)
    λ, μ, _, H = props[1], props[2], props[3], props[4]

    Fp_old  = reshape(state_old[1:9], 3, 3)
    eqps    = state_old[10]
    Fp_new  = reshape(state_new[1:9], 3, 3)
    Δεᵖ     = state_new[10] - eqps

    Fm       = Matrix{Float64}(F_sa)
    Fp_old_i = inv(Fp_old)
    Fp_new_i = inv(Fp_new)
    Fe_tr    = Fm * Fp_old_i
    Fe_new   = Fm * Fp_new_i
    Fe_new_iT = inv(Fe_new)'
    F_iT      = inv(Fm)'

    # Spectral decomposition of symmetric trial Ce
    Ce_tr = Symmetric(Fe_tr' * Fe_tr)
    eig   = eigen(Ce_tr)
    c     = eig.values      # principal squared stretches  (sorted ascending)
    Q     = eig.vectors     # columns = eigenvectors

    # Trial principal Hencky strains and Mandel stresses
    ε_tr  = 0.5 .* log.(c)
    tr_ε  = sum(ε_tr)
    m_tr  = [λ * tr_ε + 2μ * ε_tr[i] for i in 1:3]
    mdev_tr = m_tr .- sum(m_tr) / 3.0
    σvm_tr  = sqrt(1.5) * norm(mdev_tr)

    f_tr = σvm_tr - (props[3] + H * eqps)

    if f_tr ≤ 0.0
        # Elastic algorithmic moduli
        m_new = m_tr
        ĉ = [λ + 2μ * Float64(i == j) for i in 1:3, j in 1:3]
    else
        # Plastic algorithmic moduli
        N̄     = 1.5 .* mdev_tr ./ σvm_tr    # principal flow values
        m_new = m_tr .- 2μ * Δεᵖ .* N̄
        β̄     = 2μ * Δεᵖ / σvm_tr
        A     = λ + μ * β̄
        B     = 2μ * (1.0 - 1.5β̄)
        Ccoef = 2μ * β̄ - 4μ^2 / (3μ + H)
        ĉ = [A + B * Float64(i == j) + Ccoef * N̄[i] * N̄[j] for i in 1:3, j in 1:3]
    end

    # Assemble ℂ_{ABCD} in principal basis
    CC = zeros(3, 3, 3, 3)

    # Material part
    for α in 1:3, β in 1:3
        coeff = ĉ[α, β] / (2c[β])
        nα = Q[:, α]; nβ = Q[:, β]
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            CC[A, B, C, D] += coeff * nα[A] * nα[B] * nβ[C] * nβ[D]
        end
    end

    # Geometric part (α ≠ β)
    tol = 1e-10
    for α in 1:3, β in 1:3
        α == β && continue
        nα = Q[:, α]; nβ = Q[:, β]
        Δc = c[α] - c[β]
        θ = if abs(Δc) > tol * (c[α] + c[β])
            (m_new[α] - m_new[β]) / Δc
        else
            # L'Hôpital limit
            (ĉ[α, α] - ĉ[β, α]) / (2c[α])
        end
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            sAB = 0.5 * (nα[A] * nβ[B] + nα[B] * nβ[A])
            sCD = 0.5 * (nα[C] * nβ[D] + nα[D] * nβ[C])
            CC[A, B, C, D] += θ * sAB * sCD
        end
    end

    # Assemble AA_{iJkL} = −F_iT[i,L]*P[k,J]
    #   + 2 Σ_{ABCD} Fe_new_iT[i,A] CC[A,B,C,D] Fp_new_iT[B,J] Fe_tr[k,C] Fp_old_i[D,L]
    Fp_new_iT = Fp_new_i'
    AA = MArray{Tuple{3,3,3,3},Float64}(undef)
    for i in 1:3, J in 1:3, k in 1:3, L in 1:3
        val = -F_iT[i, L] * P[k, J]
        for A in 1:3, B in 1:3, C in 1:3, D in 1:3
            val += 2 * Fe_new_iT[i, A] * CC[A, B, C, D] *
                   Fp_new_iT[B, J] * Fe_tr[k, C] * Fp_old_i[D, L]
        end
        AA[i, J, k, L] = val
    end
    return SArray{Tuple{3,3,3,3},Float64,4,81}(AA)
end

# ---------------------------------------------------------------------------
# CM public API
# ---------------------------------------------------------------------------

function helmholtz_free_energy(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    F_sa = SMatrix{3,3,Float64,9}(Tensor{2,3,Float64}(∇u) + one(Tensor{2,3,Float64}))
    W, _, state_new_vec = _j2_stress(props, F_sa, Z_old)
    Z_new .= state_new_vec
    return W
end

function pk1_stress(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    F_ten = Tensor{2,3,Float64}(∇u) + one(Tensor{2,3,Float64})
    F_sa  = SMatrix{3,3,Float64,9}(F_ten)
    _, P_sa, state_new_vec = _j2_stress(props, F_sa, Z_old)
    Z_new .= state_new_vec
    return Tensor{2,3,Float64}(P_sa)
end

function material_tangent(
    ::FiniteDefJ2Plasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    F_ten = Tensor{2,3,Float64}(∇u) + one(Tensor{2,3,Float64})
    F_sa  = SMatrix{3,3,Float64,9}(F_ten)
    _, P_sa, state_new_vec = _j2_stress(props, F_sa, Z_old)
    Z_new .= state_new_vec
    AA_sa = _j2_tangent_analytical(props, F_sa, Z_old, P_sa, state_new_vec)
    # Convert SArray{Tuple{3,3,3,3}} → Tensor{4,3}
    data = ntuple(Val(81)) do lin
        l, rem = divrem(lin - 1, 27)
        k, rem = divrem(rem,      9)
        j, i   = divrem(rem,      3)
        AA_sa[i+1, j+1, k+1, l+1]
    end
    return Tensor{4,3,Float64,81}(data)
end
