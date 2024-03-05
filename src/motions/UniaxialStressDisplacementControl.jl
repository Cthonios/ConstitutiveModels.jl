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
  model, model_inputs,
  u, p
)
  F = form_deformation_gradient(motion, p[1], u)
  props, Δt, θ, state_old = model_inputs
  P, _ = pk1_stress(model, props, Δt, F, θ, state_old)

  # specialization for linear models otherwise we get singular exceptions
  if typeof(model) <: LinearElastic || 
     typeof(model) <: LinearElastoPlasticity

    Pvec = @views tovoigt(SVector, P)[2:6]
  else
    Pvec = @views tovoigt(SVector, P)[2:end]
  end
  return Pvec
end

function deformation_gradient(
  motion::Type{UniaxialStressDisplacementControl}, λ::T,
  model, model_inputs,
) where T
  # system seems singular for identity so this is a catch for now
  if λ ≈ 1.
    @show "returning 1"
    return one(Tensor{2, 3, T, 9})
  end

  x0 = SVector{8, T}([1., 1., 0., 0., 0., 0., 0., 0.])
  p  = SVector{1, Float64}((λ,))  

  f = x -> motion_objective(motion, model, model_inputs, x, p)
  x = solve_motion(f, x0)
  F = form_deformation_gradient(motion, λ, x)
  return F
end