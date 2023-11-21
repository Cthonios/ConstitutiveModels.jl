"""
"""
abstract type AbstractMotion end
"""
"""
abstract type SimpleMotion <: AbstractMotion end
"""
Returns the deformation gradient for a given motion
"""
function deformation_gradient end
function motion_objective end

@doc raw"""
Provides an analytic motion for uniaxial stress
assuming perfect incompressibility.

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda & 0 &                        0 \\
0       & \frac{1}{\sqrt{\lambda}} & 0 \\
0       & 0                        & \frac{1}{\sqrt{\lambda}}
\end{bmatrix}
```
"""
struct IsochoricUniaxialStress <: SimpleMotion
end

"""
"""
deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))

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
struct SimpleShear <: SimpleMotion
end

"""
"""
deformation_gradient(::Type{SimpleShear}, γ::T) where T <: Number = 
Tensor{2, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))

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
struct UniaxialStrain <: SimpleMotion
end

"""
"""
deformation_gradient(::Type{UniaxialStrain}, λ::T) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 

struct UniaxialStressDisplacementControl <: SimpleMotion
end

function form_deformation_gradient(
  ::Type{UniaxialStressDisplacementControl},
  λ::T, x::V
) where {T <: Number, V <: AbstractArray}
  F = Tensor{2, 3, eltype(x), 9}((
    λ,    x[8], x[7],
    x[5], x[1], x[6],
    x[4], x[3], x[2]
  ))
  return F
end

function motion_objective(
  motion::Type{UniaxialStressDisplacementControl}, 
  model::Mod, props::Props,
  u, p
) where {Mod <: ConstitutiveModel, Props <: AbstractArray}
  F    = form_deformation_gradient(motion, p[1], u)
  P    = pk1_stress(model, props, F)
  Pvec = @views tovoigt(SVector, P)[2:end]
  return Pvec
end

"""
"""
function deformation_gradient(
  motion::Type{UniaxialStressDisplacementControl},
  model::Mod, props::Props, λ::T
) where {T <: Number, Mod <: ConstitutiveModel, Props <: AbstractArray{T, 1}}

  x0      = SVector{8, T}([1., 1., 0., 0., 0., 0., 0., 0.])
  f       = (u, p) -> motion_objective(motion, model, props, u, p)
  problem = NonlinearProblem(f, x0, SVector{1, Float64}([λ]))
  sol     = solve(problem, NewtonRaphson())
  F       = form_deformation_gradient(motion, λ, sol.u)
  return F
end
