# NeoHookean Model
struct NeoHookean <: HyperelasticModel
end

function NeoHookean(inputs::Dict{String, F}) where F <: AbstractFloat
  @assert "bulk modulus" in keys(inputs)
  @assert "shear modulus" in keys(inputs)
  @assert inputs["shear modulus"] > 0.0

  props = F[
    inputs["bulk modulus"], 
    inputs["shear modulus"]
  ]
  return NeoHookean(), props
end
 
function strain_energy_density(::NeoHookean, F::T, props::V) where {T <: Tensor{2, 3, <:Number}, V <: AbstractArray{<:Number, 1}}
  K, G    = props[1], props[2]
  J       = det(F)
  I_1_bar = tr(J^(-2. / 3.) * tdot(F))
  W_vol   = 0.5 * K * (0.5 * (J^2 - 1) - log(J))
  W_dev   = 0.5 * G * (I_1_bar - 3.)
  return W_vol + W_dev
end

function strain_energy_density(::NeoHookean, C::T, props::V) where {T <: SymmetricTensor{2, 3, <:Number}, V <: AbstractArray{<:Number, 1}}
  K, G    = props[1], props[2]
  I_3     = det(C)
  I_1_bar = tr(I_3^(-1. / 3.) * C)
  W_vol   = 0.5 * K * (0.5 * (I_3 - 1) - 0.5 * log(I_3))
  W_dev   = 0.5 * G * (I_1_bar - 3.)
  return W_vol + W_dev
end