abstract type AbstractSimpleMotion end

struct MotionNewtonSolverSettings
  max_iter::Int32
  abs_tol::Float64
  rel_tol::Float64
end

function MotionNewtonSolverSettings()
  return MotionNewtonSolverSettings(Int32(20), 1e-12, 1e-12)
end

function solve(f, x0, model)
  settings = MotionNewtonSolverSettings()
  R0 = 1.0e6
  x = x0
  for n in 1:settings.max_iter
    R = f(x, model)

    if n == 1
      R0 = norm(R)
    end

    if norm(R) < settings.abs_tol
      break
    end

    if norm(R) / R0 < settings.rel_tol
      break
    end

    K = FiniteDiff.finite_difference_jacobian(z -> f(z, model), x)
    # display(R)
    # display(K)
    Δx = -K \ R
    x = x + Δx

    if norm(Δx) < settings.abs_tol
      break
    end

    if n == settings.max_iter
      @warn "Reached maximum Newton iterations. Be careful"
    end
  end

  return x
end

# stuff for nonlinear motions
function motion_objective(motion, u, p, model, props, args...)
  ∇u = _displacement_gradient(motion, p[1], u)
  P, _ = pk1_stress(model, props, ∇u, args...)
  return _extract_stress(motion, model, P)
end

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
struct IsochoricUniaxialStress <: AbstractSimpleMotion
end

displacement_gradient(::Type{IsochoricUniaxialStress}, λ::T) where T <: Number = 
Tensor{2, 3, T, 9}((λ - 1, 0., 0., 0., 1. / sqrt(λ) - 1, 0., 0., 0., 1. / sqrt(λ) - 1))

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
struct SimpleShear <: AbstractSimpleMotion
end

displacement_gradient(::Type{SimpleShear}, γ::T) where T <: Number = 
Tensor{2, 3, T, 9}((0., 0., 0., γ, 0., 0., 0., 0., 0.))

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
struct UniaxialStrain <: AbstractSimpleMotion
end

displacement_gradient(::Type{UniaxialStrain}, λ::T, args...) where T <: Number = 
Tensor{2, 3, T, 9}((λ - 1, 0., 0., 0., 0., 0., 0., 0., 0.))

# Stress control
struct UniaxialStressDisplacementControl <: AbstractSimpleMotion
end

function _displacement_gradient(::Type{UniaxialStressDisplacementControl}, λ::T, x) where T <: Number
  ∇u = Tensor{2, 3, eltype(x), 9}((
    λ - 1, x[8], x[7],
    x[5],  x[1], x[6],
    x[4],  x[3], x[2]
  ))
  return ∇u
end

function _extract_stress(::Type{UniaxialStressDisplacementControl}, model, P)
  if (isa(model, Hyperelastic) && isa(model.strain_energy, LinearElastic)) ||
     isa(model, LinearElastoPlastic)
    return tovoigt(SVector, P)[2:6]
  else
    return tovoigt(SVector, P)[2:end]
  end
end

function displacement_gradient(motion::Type{UniaxialStressDisplacementControl}, λ::T, args...) where T <: Number
  x0 = SVector{8, T}([0., 0., 0., 0., 0., 0., 0., 0.])
  p = SVector{1, Float64}((λ,))
  f = (x, m) -> motion_objective(motion, x, p, m, args[2:end]...)
  x = solve(f, x0, args[1])
  return _displacement_gradient(motion, λ, x)
end


# Helpers
struct MaterialState{T1, T2, T3, T4, T5, T6}
  F::T1
  ψ::T2
  P::T3
  σ::T4
  A::T5
  state::T6
end

function MaterialState(model, inputs, motion, ps...)
  model = ConstitutiveModel(model, inputs)
  props = initialize_properties(model, inputs)
  state_old = initialize_state(model)
  mat_states = Vector{MaterialState}(undef, 0)

  # TODO hook this up to constructor eventually
  θ = 0.0
  Δt = 1.0

  for p in zip(ps...)
    ∇u = displacement_gradient(motion, p..., model, props, θ, state_old, Δt)
    F = deformation_gradient(∇u)
    ψ, state = helmholtz_free_energy(model, props, ∇u, θ, state_old, Δt)
    P, state = pk1_stress(model, props, ∇u, θ, state_old, Δt)
    σ, state = cauchy_stress(model, props, ∇u, θ, state_old, Δt)
    A, state = material_tangent(model, props, ∇u, θ, state_old, Δt)
    push!(mat_states, MaterialState(F, ψ, P, σ, A, state))
    state_old = state
  end

  return mat_states
end
