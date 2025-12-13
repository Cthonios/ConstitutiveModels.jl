module Tensors

using DifferentiationInterface
using Enzyme
using LinearAlgebra
using StaticArrays

const Tensor1{T} = SVector{3, T} where T <: Number
const Tensor2{T} = SMatrix{3, 3, T, 9} where T <: Number
const Tensor3{T} = SArray{Tuple{3, 3, 3}, T, 3, 27} where T <: Number
const Tensor4{T} = SArray{Tuple{3, 3, 3, 3}, T, 4, 81} where T <: Number

function tdot(F::Tensor2)
    return F' * F
end

struct NeoHookean
end

function energy(::NeoHookean, props, Δt, ∇u, θ, Z_old, Z_new)
    # unpack properties
    κ, μ = props[1], props[2]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
    J       = det(F)

    J_m_13  = 1. / cbrt(J)
    J_m_23  = J_m_13 * J_m_13
    I_1_bar = tr(J_m_23 * tdot(F))
    ψ_vol = 0.5 * κ * (0.5 * (J * J - one(eltype(J))) - log(J))
    ψ_dev = 0.5 * μ * (I_1_bar - 3.)
    ψ     = ψ_vol + ψ_dev
    return ψ
end

function my_func(model, props, Δt, ∇u, θ, Z_old, Z_new)
    return tr(∇u)
end

# function denerg

function test()
    model = NeoHookean()
    props = SVector{2, Float64}(100., 1.)
    Δt = 0.
    ∇u = one(Tensor2)
    # d∇u = zero(Tensor2)
    d∇u = Tensor2(ones(9))
    θ = 0.
    Z_old = zeros(0)
    dZ_old = zeros(0)
    Z_new = zeros(0)
    dZ_new = zeros(0)

    @info "Energy calculation"
    @time energy(model, props, Δt, ∇u, θ, Z_old, Z_new)

    @time autodiff(
        Reverse,
        energy,
        Const(model),
        Const(props),
        Const(Δt),
        Active(∇u),
        Const(θ),
        Duplicated(Z_old, dZ_old),
        Duplicated(Z_new, dZ_new)
    )

    # @time autodiff(
    #     Reverse,
    #     my_func,
    #     Active(∇u)
    # )
    # @time autodiff(
    #     Reverse,
    #     my_func,
    #     Const(model),
    #     Const(props),
    #     Const(Δt),
    #     Active(∇u),
    #     Const(θ),
    #     Duplicated(Z_old, dZ_old),
    #     Duplicated(Z_new, dZ_new)
    # )
    @time out = autodiff(
        Forward,
        energy,
        Duplicated,
        Const(model),
        Const(props),
        Const(Δt),
        # Active(∇u),
        # ∇u,
        Duplicated(∇u, d∇u),
        Const(θ),
        Duplicated(Z_old, dZ_old),
        Duplicated(Z_new, dZ_new)
    )
    d∇u
    out

    # backend = AutoEnzyme(mode=Enzyme.Reverse)
    # func = z -> energy(model, props, Δt, z, θ, Z_old, Z_new)
    # prep = prepare_gradient(func, backend, ∇u)
    # @time df = DifferentiationInterface.gradient(func, prep, backend, ∇u)
    # @time df = DifferentiationInterface.gradient(func, prep, backend, ∇u)

end

end # module
