module ConstitutiveModelsEnzymeExt

using ConstitutiveModels
using Enzyme
using StaticArrays
using Tensors

"""
useful for hyperelastic models
"""
function ConstitutiveModels.helmholtz_free_energy_sensitivies(
  ::Enzyme.ReverseMode{false},
  model::Mod, props::V1, F::M
) where {
  Mod <: MechanicalModel,
  V1  <: AbstractArray{<:Number, 1},
  M   <: AbstractArray{<:Number, 2}
}

  grads = autodiff(Reverse, helmholtz_free_energy, model, Active(props), Active(F))
  return grads[1][2], grads[1][3]
end

function ConstitutiveModels.helmholtz_free_energy_sensitivies(
  ::Enzyme.ReverseMode{false},
  model::Mod, props::V1, F::Tensor{2, 3, T, 9}, state_old::V2
) where {
  Mod <: MechanicalModel,
  V1  <: AbstractArray{<:Number, 1},
  T   <: Number,
  V2  <: AbstractArray{<:Number, 1}
}

  grads = autodiff(Reverse, helmholtz_free_energy, model, Active(props), Active(F), Active(state_old))
  return grads[1][2], grads[1][3], grads[1][4]
end

function ConstitutiveModels.helmholtz_free_energy_sensitivies(
  ::Enzyme.ReverseMode{false},
  model::Mod, props::V1, ε::SymmetricTensor{2, 3, T, 6}, state_old::V2
) where {
  Mod <: MechanicalModel,
  V1  <: AbstractArray{<:Number, 1},
  T   <: Number,
  V2  <: AbstractArray{<:Number, 1}
}

  grads = autodiff(Reverse, helmholtz_free_energy, model, Active(props), Active(ε), Active(state_old))
  return grads[1][2], grads[1][3], grads[1][4]
end

end # module