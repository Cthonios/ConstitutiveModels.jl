
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

# function form_deformation_gradient(
#   ::Type{UniaxialStressDisplacementControl},
#   λ::T, x::V, type::Type{MMatrix}
# ) where {T <: Number, V <: AbstractArray}
#   F = MMatrix{3, 3, eltype(x), 9}((
#     λ,    x[8], x[7],
#     x[5], x[1], x[6],
#     x[4], x[3], x[2]
#   ))
#   return F
# end

# function form_deformation_gradient(
#   ::Type{UniaxialStressDisplacementControl},
#   λ::T, x::V, type::Type{SMatrix}
# ) where {T <: Number, V <: AbstractArray}
#   F = SMatrix{3, 3, eltype(x), 9}((
#     λ,    x[8], x[7],
#     x[5], x[1], x[6],
#     x[4], x[3], x[2]
#   ))
#   return F
# end

# function form_deformation_gradient(
#   ::Type{UniaxialStressDisplacementControl},
#   λ::T, x::V, type::Type{Matrix}
# ) where {T <: Number, V <: AbstractArray}
#   F = T[
#     λ    x[8] x[7];
#     x[5] x[1] x[6];
#     x[4] x[3] x[2]
#   ]
#   return F
# end

# TODO need to specialize to regular ole matrices
# with form_deformation_gradient

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

function motion_objective(
  motion::Type{UniaxialStressDisplacementControl}, 
  model::Mod, props::Props, state_old,
  u, p
) where {Mod <: ConstitutiveModel, Props <: AbstractArray}
  F = form_deformation_gradient(motion, p[1], u)
  # P    = pk1_stress(model, props, F)
  # P    = pk1_stress!(model, props, F, state_old, state_new)
  P, props, state_new = pk1_stress(model, props, F, state_old)
  Pvec = @views tovoigt(SVector, P)[2:end]
  return Pvec

end

# TODO add solver options
"""
"""
function deformation_gradient(
  motion::Type{UniaxialStressDisplacementControl},
  model::Mod, props::Props, λ::T, type::Type
) where {T <: Number, Mod <: ConstitutiveModel, Props <: AbstractArray{T, 1}}

  # system seems singular for identity so this is a catch for now
  if λ ≈ 1.
    return one(Tensor{2, 3, T, 9})
  end

  x0 = SVector{8, T}([1., 1., 0., 0., 0., 0., 0., 0.])
  p  = SVector{1, Float64}((λ,))  
  f  = x -> motion_objective(motion, model, props, x, p)
  x  = solve_motion(f, x0)
  F  = form_deformation_gradient(motion, λ, x)
  return F
end

function deformation_gradient(
  motion::Type{UniaxialStressDisplacementControl},
  model::Mod, props::Props, state_old,
  λ::T, type::Type
) where {T <: Number, Mod <: ConstitutiveModel, Props <: AbstractArray{T, 1}}

  # system seems singular for identity so this is a catch for now
  if λ ≈ 1.
    return one(Tensor{2, 3, T, 9})
  end

  x0 = SVector{8, T}([1., 1., 0., 0., 0., 0., 0., 0.])
  p  = SVector{1, Float64}((λ,))  
  # f  = x -> motion_objective(motion, model, props, x, p)
  # f  = (x, s) -> motion_objective!(motion, model, props, state_old, s, x, p)
  # f  = x -> motion_objective!(motion, model, props, state_old, state_new, x, p)
  f  = x -> motion_objective(motion, model, props, state_old, x, p)
  # x  = solve_motion(f, x0, state_new)
  x  = solve_motion(f, x0)
  F  = form_deformation_gradient(motion, λ, x)
  return F
end