module ConstitutiveModelsEnzymeExt

using ConstitutiveModels
using Enzyme
using StaticArrays
using Tensors

# general AD stuff
function ConstitutiveModels.pk1_stress(
  ::Enzyme.ReverseMode{false}, 
  model::Mod, props::V1, F::M, state::V2
) where {Mod <: HyperelasticModel, 
         V1  <: AbstractArray{<:Number, 1}, 
         M   <: AbstractArray{<:Number, 2},
         V2  <: AbstractArray{<:Number, 1}}
  
  # grads = autodiff(Reverse, energy!, Active, model, props, Active(F), state)
  # return grads
  grads = autodiff(Reverse, energy, Active, model, props, Active(F), state)
  return grads[1][3]
end

function ConstitutiveModels.cauchy_stress(
  ::Enzyme.ReverseMode{false}, 
  model::Mod, props::V1, F::M, state::V2
) where {Mod <: HyperelasticModel, 
         V1  <: AbstractArray{<:Number, 1}, 
         M   <: AbstractArray{<:Number, 2},
         V2  <: AbstractArray{<:Number, 1}}

  grads = autodiff(Reverse, energy, Active, model, props, Active(F), state)
  P = grads[1][3]
  J = det(F)
  σ = J^-1 * dot(P, F') |> symmetric
end

end # module