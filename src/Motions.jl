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
deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T, type::Type{<:Tensor}) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))
deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T, type::Type{<:MMatrix}) where T <: Number = 
MMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))
deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T, type::Type{<:SMatrix}) where T <: Number = 
SMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))
deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T, type::Type{<:Matrix}) where T <: Number = 
T[λ 0. 0.; 0. 1. / sqrt(λ) 0.; 0. 0. 1. / sqrt(λ)]

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
deformation_gradient(::Type{SimpleShear}, γ::T, ::Type{<:Tensor}) where T <: Number = 
Tensor{2, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))
deformation_gradient(::Type{SimpleShear}, γ::T, ::Type{<:MMatrix}) where T <: Number = 
MMatrix{3, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))
deformation_gradient(::Type{SimpleShear}, γ::T, ::Type{<:SMatrix}) where T <: Number = 
SMatrix{3, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))
deformation_gradient(::Type{SimpleShear}, γ::T, type::Type{<:Matrix}) where T <: Number = 
T[1. γ 0.; 0. 1. 0.; 0. 0. 1.]

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
deformation_gradient(::Type{UniaxialStrain}, λ::T, ::Type{<:Tensor}) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 
deformation_gradient(::Type{UniaxialStrain}, λ::T, ::Type{<:MMatrix}) where T <: Number = 
MMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 
deformation_gradient(::Type{UniaxialStrain}, λ::T, ::Type{<:SMatrix}) where T <: Number = 
SMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 
deformation_gradient(::Type{UniaxialStrain}, λ::T, type::Type{<:Matrix}) where T <: Number = 
T[λ 0. 0.; 0. 1. 0.; 0. 0. 1.]

struct UniaxialStressDisplacementControl <: SimpleMotion
end

function form_deformation_gradient(
  ::Type{UniaxialStressDisplacementControl},
  λ::T, x::V, ::Type{<:Tensor}
) where {T <: Number, V <: AbstractArray}
  F = Tensor{2, 3, eltype(x), 9}((
    λ,    x[8], x[7],
    x[5], x[1], x[6],
    x[4], x[3], x[2]
  ))
  return F
end

function form_deformation_gradient(
  ::Type{UniaxialStressDisplacementControl},
  λ::T, x::V, type::Type{MMatrix}
) where {T <: Number, V <: AbstractArray}
  F = MMatrix{3, 3, eltype(x), 9}((
    λ,    x[8], x[7],
    x[5], x[1], x[6],
    x[4], x[3], x[2]
  ))
  return F
end

function form_deformation_gradient(
  ::Type{UniaxialStressDisplacementControl},
  λ::T, x::V, type::Type{SMatrix}
) where {T <: Number, V <: AbstractArray}
  F = SMatrix{3, 3, eltype(x), 9}((
    λ,    x[8], x[7],
    x[5], x[1], x[6],
    x[4], x[3], x[2]
  ))
  return F
end

function form_deformation_gradient(
  ::Type{UniaxialStressDisplacementControl},
  λ::T, x::V, type::Type{Matrix}
) where {T <: Number, V <: AbstractArray}
  F = T[
    λ    x[8] x[7];
    x[5] x[1] x[6];
    x[4] x[3] x[2]
  ]
  return F
end

# TODO need to specialize to regular ole matrices
# with form_deformation_gradient

function motion_objective(
  motion::Type{UniaxialStressDisplacementControl}, 
  model::Mod, props::Props, state, type::Type,
  u, p
) where {Mod <: ConstitutiveModel, Props <: AbstractArray}
  F        = form_deformation_gradient(motion, p[1], u, type)
  P, state = pk1_stress(model, props, F, state)
  Pvec     = @views tovoigt(SVector, P)[2:end]
  return Pvec
end

"""
"""
function deformation_gradient(
  motion::Type{UniaxialStressDisplacementControl},
  model::Mod, props::Props, state, λ::T, type::Type
) where {T <: Number, Mod <: ConstitutiveModel, Props <: AbstractArray{T, 1}}

  # system seems singular for identity so this is a catch for now
  if λ ≈ 1.
    if type <: Tensor
      return one(SymmetricTensor{2, 3, T, 6})
    elseif type <: MMatrix
      return MMatrix{3, 3, T, 9}((1., 0., 0., 0., 1., 0., 0., 0., 1.))
    elseif type <: SMatrix
      return SMatrix{3, 3, T, 9}((1., 0., 0., 0., 1., 0., 0., 0., 1.))
    elseif type <: Matrix
      return T[1. 0. 0.; 0. 1. 0.; 0. 0. 1.]
    end
  end

  x0      = SVector{8, T}([1., 1., 0., 0., 0., 0., 0., 0.])
  f       = (u, p) -> motion_objective(motion, model, props, state, Tensor, u, p)
  problem = NonlinearProblem(f, x0, SVector{1, Float64}([λ]))
  sol     = solve(problem, NewtonRaphson())
  F       = form_deformation_gradient(motion, λ, sol.u, type)
  return F
end


# default 
const DefaultMotionType = MMatrix
deformation_gradient(motion::Type{<:SimpleMotion}, λ::T; type::Type = DefaultMotionType) where T <: Number = 
deformation_gradient(motion, λ, type)
deformation_gradient(motion::Type{<:SimpleMotion}, model::Mod, props::Props, state, λ::T; type::Type = DefaultMotionType) where {Mod <: MechanicalModel, Props <: AbstractArray, T <: Number} =
deformation_gradient(motion, model, props, state, λ, type)
