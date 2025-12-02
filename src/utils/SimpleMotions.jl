"""
$(TYPEDEF)
"""
abstract type AbstractSimpleMotion{T <: Number} end
Base.eltype(::AbstractSimpleMotion{T}) where T = T

"""
$(TYPEDSIGNATURES)
"""
function displacement_gradient(motion::AbstractSimpleMotion{T}, args...) where T
    F = deformation_gradient(motion, args...)
    I = one(Tensor{2, 3, eltype(motion), 9})
    return F - I
end

@doc raw"""
Provides an analytic motion for uniaxial strain

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda_1 & 0         & 0 \\
0         & \lambda_2 & 0 \\
0         & 0         & 1
\end{bmatrix}
```
"""
struct BiaxialStrain{T} <: AbstractSimpleMotion{T}
end

function BiaxialStrain()
    return BiaxialStrain{Float64}()
end

function deformation_gradient(motion::BiaxialStrain{T}, λ1::T, λ2::T) where T
    Tensor{2, 3, eltype(motion), 9}((
        λ1, 0., 0.,
        0., λ2, 0.,
        0., 0., 1.
    ))
end

@doc raw"""
Provides an analytic motion for uniaxial strain

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda & 0                 & 0 \\
0       & \frac{1}{\lambda} & 0 \\
0       & 0                 & \frac{1}{\lambda}
\end{bmatrix}
```
"""
struct IsochoricUniaxialStress{T} <: AbstractSimpleMotion{T}
end

function IsochoricUniaxialStress()
    return IsochoricUniaxialStress{Float64}()
end

function deformation_gradient(motion::BiaxialStrain{T}, λ::T) where T
    Tensor{2, 3, eltype(motion), 9}((
        λ,  0.,           0.,
        0., 1. / sqrt(λ), 0.,
        0., 0.,           1. / sqrt(λ)
    ))
end

@doc raw"""
Provides an analytic motion for uniaxial strain

This is
```math
\mathbf{F} = \frac{1}{2}\begin{bmatrix}
\left(\lambda + \lambda^{-1}\right) & \left(\lambda - \lambda^{-1}\right) & 0 \\
\left(\lambda - \lambda^{-1}\right) & \left(\lambda + \lambda^{-1}\right) & 0 \\
0                                   & 0                                   & 2
\end{bmatrix}
```
"""
struct PureShearStrain{T} <: AbstractSimpleMotion{T}
end

function PureShearStrain()
    return PureShearStrain{Float64}()
end

function deformation_gradient(motion::PureShearStrain{T}, λ::T) where T
    0.5 * Tensor{2, 3, eltype(motion), 9}((
        λ + 1. / λ, λ - 1. / λ, 0.,
        λ - 1. / λ, λ + 1. / λ, 0.,
        0.,         0.,         2.
    ))
end

@doc raw"""
Provides an analytic motion for simple shear.

This is
```math
\mathbf{F} = \begin{bmatrix}
1 & \gamma & 0 \\
0 & 1      & 0 \\
0 & 0      & 1
\end{bmatrix}
```
"""
struct SimpleShear{T} <: AbstractSimpleMotion{T}
end

function SimpleShear()
    return SimpleShear{Float64}()
end

function deformation_gradient(motion::SimpleShear{T}, γ::T) where T
    Tensor{2, 3, eltype(motion), 9}((
        1., 0., 0.,
        γ,  1., 0.,
        0., 0., 1.
    ))
end

@doc raw"""
Provides an analytic motion for uniaxial strain

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda & 0 & 0 \\
0       & 1 & 0 \\
0       & 0 & 1
\end{bmatrix}
```
"""
struct UniaxialStrain{T} <: AbstractSimpleMotion{T}
end

function UniaxialStrain()
    return UniaxialStrain{Float64}()
end

function deformation_gradient(motion::UniaxialStrain{T}, λ::T) where T
    Tensor{2, 3, eltype(motion), 9}((
        λ,  0., 0.,
        0., 1., 0.,
        0., 0., 1.
    ))
end

"""
Document me
"""
struct UniaxialStressDisplacementControl{T} <: AbstractSimpleMotion{T}
end

function UniaxialStressDisplacementControl()
    return UniaxialStressDisplacementControl{Float64}()
end

function _motion_objective(
    ::UniaxialStressDisplacementControl{T},
    x, p
) where T <: Number
    λ, model, model_inputs = p
    ∇u = Tensor{2, 3, eltype(x), 9}((
        λ - 1., 0., 0.,
        0.,     x,  0.,
        0.,     0., x
    ))
    props, Δt, θ, Z_old, Z_new = model_inputs
    σ = pk1_stress(model, props, Δt, ∇u, θ, Z_old, Z_new)
    return σ[2, 2]
end

function deformation_gradient(
    motion::UniaxialStressDisplacementControl{T},
    λ::T,
    model,
    model_inputs...
) where T <: Number
    ∇u = displacement_gradient(motion, λ, model, model_inputs...)
    F = ∇u + one(∇u)
    return F
end

function displacement_gradient(
    motion::UniaxialStressDisplacementControl{T},
    λ::T,
    model,
    model_inputs...
) where T <: Number

    if λ ≈ 1
        return zero(Tensor{2, 3, T, 9})
    end

    p  = (λ, model, model_inputs...)
    f = x -> _motion_objective(motion, x, p)
    x = find_zero(f, 1.)
    ∇u = Tensor{2, 3, eltype(x), 9}((
        λ - 1., 0., 0.,
        0.,     x,  0.,
        0.,     0., x
    ))
    return ∇u
end

function simulate_material_point(
    f,
    model::AbstractConstitutiveModel{NP, NS},
    props, Δt,
    θ, Z_old, Z_new,
    motion::AbstractSimpleMotion,
    args...
) where {NP, NS}
    

    ∇u = zero(Tensor{2, 3, typeof(θ), 9})
    result_type = typeof(f(model, props, Δt, ∇u, θ, Z_old, Z_new))

    ∇us = Tensor{2, 3, typeof(θ), 9}[]
    results = result_type[]
    state_news = typeof(Z_new)[]

    for params in zip(args...)
        if isa(motion, UniaxialStressDisplacementControl)
            # @assert false
            ∇u = displacement_gradient(
                motion, params..., 
                model, (props, Δt, θ, Z_old, Z_new)
            )
        else
            ∇u = displacement_gradient(motion, params...)
        end
        a = f(model, props, Δt, ∇u, θ, Z_old, Z_new)
        push!(∇us, ∇u)
        push!(results, a)
        push!(state_news, copy(Z_new))
    end

    # as = map(x -> x[1], results)
    # Zs = map(x -> x[2], results)

    return ∇us, results, state_news
    # return ∇us, as, Zs
end