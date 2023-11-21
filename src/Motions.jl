abstract type AbstractMotion end
abstract type SimpleMotion <: AbstractMotion end

function deformation_gradient end
function motion_objective end

struct IsochoricUniaxialStress <: SimpleMotion
end

deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))

struct SimpleShear <: SimpleMotion
end

deformation_gradient(::Type{SimpleShear}, γ::T) where T <: Number = 
Tensor{2, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))

struct UniaxialStrain <: SimpleMotion
end

deformation_gradient(::Type{UniaxialStrain}, λ::T) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 

struct UniaxialStressDisplacementControl <: SimpleMotion
end

function form_deformation_gradient(
  ::Type{UniaxialStressDisplacementControl},
  λ::T, x::SVector{8, T}
) where T <: Number
  F = Tensor{2, 3, Float64, 9}((
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

function deformation_gradient(
  motion::Type{UniaxialStressDisplacementControl},
  model::Mod, props::Props, λ::T
) where {T <: Number, Mod <: ConstitutiveModel, Props <: AbstractArray{T, 1}}

  x0      = SVector{8, T}([1., 1., 0., 0., 0., 0., 0., 0.])
  f       = (u, p) -> motion_objective(motion, model, props, u, p)
  problem = NonlinearProblem(f, x0, SVector{1, Float64}([λ]))
  sol     = solve(problem, SimpleNonlinearSolve.SimpleNewtonRaphson(; autodiff=false))
  F       = form_deformation_gradient(motion, λ, sol.u)
  return F
end
