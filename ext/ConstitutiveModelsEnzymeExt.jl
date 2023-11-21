module ConstitutiveModelsEnzymeExt

using ConstitutiveModels
using Enzyme
using Tensors

# general AD stuff
function ConstitutiveModels.pk1_stress(
  ::Enzyme.ReverseMode{false}, 
  model::M, props::V, F::T
) where {M <: HyperelasticModel, 
         T <: Union{Tensor{2, 3, <:Number, 9}, SymmetricTensor{2, 3, <:Number, 6}}, 
         V <: AbstractArray{<:Number, 1}}
  
  # annoying for right now..., maybe write a better wrapper
  autodiff(Reverse, energy, Active, model, props, Active(F))[1][3]
end

function ConstitutiveModels.pk2_stress(
  ::Enzyme.ReverseMode{false}, 
  model::M, props::V, C::T
) where {M <: HyperelasticModel, 
         T <: SymmetricTensor{2, 3, <:Number, 6}, 
         V <: AbstractArray{<:Number, 1}}
  
  # annoying for right now..., maybe write a better wrapper
  2. * autodiff(Reverse, energy, Active, model, props, Active(C))[1][3]
end

function ConstitutiveModels.cauchy_stress(
  ::Enzyme.ReverseMode{false}, 
  model::M, props::V, F::T
) where {M <: HyperelasticModel, 
         T <: Union{Tensor{2, 3, <:Number, 9}, SymmetricTensor{2, 3, <:Number, 6}}, 
         V <: AbstractArray{<:Number, 1}}
  
  # annoying for right now..., maybe write a better wrapper
  if typeof(model) <: ConstitutiveModels.LinearElastic
    I = one(SymmetricTensor{2, 3, eltype(F), 6})
    ∇u = F - I
    ε = 0.5 * (∇u + ∇u')
    σ = autodiff(Reverse, ConstitutiveModels.energy_enzyme_specialization, Active, model, props, Active(ε))[1][3]
    return σ
  else
    J = det(F)
    P = autodiff(Reverse, energy, Active, model, props, Active(F))[1][3]
    return J^-1 * symmetric(dot(P, F'))
  end
end

# function ConstitutiveModels.cauchy_stress(
#   ::Enzyme.ReverseMode{false},
#   model::ConstitutiveModels.LinearElastic, props::V, F::Tensor{2, 3, <:Number, 9}
# ) where {V <: AbstractArray}

#   @show "here"
#   I = one(SymmetricTensor{2, 3, eltype(F), 6})
#   ∇u = F - I
#   ε  = 0.5 * symmetric(∇u + ∇u')
#   σ = autodiff(Reverse, energy, ACtive, model, props, Active(ε))[1][3]
#   return σ
# end

# function ConstitutiveModels.strain_energy_density_and_property_adjoints(
#   model::M, F::T, props::V
# ) where {M <: HyperelasticModel, T <: Tensor{2, 3, <:Number, 9}, V <: AbstractArray{<:Number, 1}}

#   grads = autodiff(ReverseWithPrimal, strain_energy_density, Active, Const(model), Const(F), Active(props))
#   return grads[2], grads[1][3]
# end

# function ConstitutiveModels.pk1_stress_and_property_adjoints(
#   model::M, F::T, props::V
# ) where {M <: HyperelasticModel, T <: Tensor{2, 3, <:Number, 9}, V <: AbstractArray{<:Number, 1}}

#   grads = autodiff(Reverse, strain_energy_density, Active, model, Active(F), Active(props))
#   return grads[1][2], grads[1][3]
# end

end # module