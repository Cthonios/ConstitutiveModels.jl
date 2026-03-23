"""
$(TYPEDEF)
"""
struct Hencky <: AbstractHyperelastic{2, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::Hencky, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    return [elastic_props.κ, elastic_props.μ]
end

"""
``
\\psi = \\frac{1}{2}\\kappa tr\\left(\\mathbf{E}\\right)^2 +
        \\mu\\mathbf{E_0}:\\mathbf{E_0}
``
$(TYPEDSIGNATURES)
"""
function strain_energy_density(::Hencky, props, F, θ)
    κ, μ = props[1], props[2]

    # kinematics
    C = tdot(F)
    E = 0.5 * log(C)
    trE = tr(E)
    devE = dev(E)

    # constitutive
    ψ = 0.5κ * trE * trE + μ * dcontract(devE, devE)
    return ψ
end

"""
$(TYPEDSIGNATURES)
"""
function pk2_stress(::Hencky, props, F, θ)
    κ, μ = props[1], props[2]

    # kinematics
    C = tdot(F)
    E = 0.5 * log(C)
    trE = tr(E)
    devE = dev(E)

    # constitutive
    dψdE = κ * trE * one(C) + 2μ * devE
    S = inv(C) ⋅ dψdE
    return S
end

"""
$(TYPEDSIGNATURES)
"""
function pk1_stress(model::Hencky, props, F, θ)
    S = pk2_stress(model, props, F, θ)
    return F ⋅ S
end

"""
$(TYPEDSIGNATURES)
"""
function material_tangent(::Hencky, props, F::Tensor{2, 3, T, 9}, θ) where T <: Number
    κ, μ = props[1], props[2]

    # kinematics
    C  = tdot(F)
    IC = inv(C)
    I2 = one(C)

    E       = 0.5 * log(C)
    trE     = tr(E)
    devE    = dev(E)
    A       = κ * trE * I2 + 2μ * devE
    S       = IC ⋅ A

    # H = ∂²ψ/∂E²
    IxI  = I2 ⊗ I2
    Isym = 0.5 * (otimesu(I2, I2) + otimesl(I2, I2))
    H    = κ * IxI + 2μ * (Isym - (1/3) * IxI)

    dist = norm(C - I2)
    dlogCdC = dist < sqrt(eps(real(T))) ? Isym : dlog(C)

    # ∂A_{ij}/∂C_{kl} = 0.5 * H_{ijmn} * dlogCdC_{mnkl}
    dAdC = Tensor{4, 3, T, 81}() do i, j, k, l
        s = zero(T)
        @inbounds for m in 1:3, n in 1:3
            s += H[i, j, m, n] * dlogCdC[m, n, k, l]
        end
        0.5 * s
    end

    # term1: IC_{im} ∂A_{mj}/∂C_{kl}
    term1 = Tensor{4, 3, T, 81}() do i, j, k, l
        s = zero(T)
        @inbounds for m in 1:3
            s += IC[i, m] * dAdC[m, j, k, l]
        end
        s
    end

    # term2: ∂(C⁻¹)_{im}/∂C_{kl} * A_{mj} = -½(IC_{ik}S_{lj} + IC_{il}S_{kj})
    term2 = Tensor{4, 3, T, 81}() do i, j, k, l
        -0.5 * (IC[i, k] * S[l, j] + IC[i, l] * S[k, j])
    end

    # ℂ = 2∂S/∂C in [i,j,k,l] = [m,J,n,L] index order.
    # _convect_tangent reads ℂ[m,j,l,n], i.e. it expects [m,J,L,n] order,
    # so we must swap indices 3 and 4 before passing in.
    dSdC = term1 + term2
    ℂ = Tensor{4, 3, T, 81}() do i, j, k, l
        2 * dSdC[i, j, l, k]   # swap k↔l to match _convect_tangent's [m,J,L,n]
    end

    return _convect_tangent(ℂ, S, F)
end

p_wave_modulus(::Hencky, props) = props[1] + 4 * props[2] / 3
