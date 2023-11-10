abstract type MechanicalModel <: ConstitutiveModel end
abstract type HyperelasticModel <: ConstitutiveModel end

# general AD stuff
pk1_stress(model::M, F::T, props::V) where {M <: HyperelasticModel, T <: Tensor{2, 3, <:Number}, V <: AbstractArray{<:Number, 1}} = 
Tensors.gradient(z -> strain_energy_density(model, z, props), F)
pk2_stress(model::M, C::T, props::V) where {M <: HyperelasticModel, T <: SymmetricTensor{2, 3, <:Number}, V <: AbstractArray{<:Number, 1}} = 
2. * Tensors.gradient(z -> strain_energy_density(model, z, props), C)

material_tangent(model::M, F::T, props::V) where {M <: HyperelasticModel, T <: Tensor{2, 3, <:Number}, V <: AbstractArray{<:Number, 1}} =
Tensors.hessian(z -> strain_energy_density(model, z, props), F)

# acoustic_tensor(model::M, F::T, props::V) where {M <: HyperelasticModel, T <: Tensor{2, 3, <:Number}, V <: AbstractArray{<:Number, 1}} =

property_adjoint(model::M, F::T, props::V) where {M <: HyperelasticModel, T <: Tensor{2, 3, <:Number}, V <: AbstractArray{<:Number, 1}} = 
ForwardDiff.gradient(z -> strain_energy_density(model, F, z), props)

# modles to include
include("mechanical_models/NeoHookean.jl")
