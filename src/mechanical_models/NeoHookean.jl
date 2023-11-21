struct NeoHookean <: HyperelasticModel{2, 0}
end

function NeoHookean(inputs::D) where D <: Dict
  @assert "bulk modulus" in keys(inputs)
  @assert "shear modulus" in keys(inputs)
  # @assert inputs["shear modulus"] > 0.0

  bulk_modulus  = inputs["bulk modulus"]
  shear_modulus = inputs["shear modulus"]

  # TODO maybe make this a helper method
  if typeof(shear_modulus) <: Quantity
    @assert ustrip(shear_modulus) > 0.0
  else
    @assert shear_modulus > 0.0
  end

  props = @SVector [
    bulk_modulus,
    shear_modulus
  ]
  return NeoHookean(), props
end
 
function energy(::NeoHookean, props::V, F::T) where {V <: AbstractArray{<:Number, 1}, T <: Tensor{2, 3, <:Number}}
  K, G    = props[1], props[2]
  J       = det(F)
  I_1_bar = tr(J^(-2. / 3.) * tdot(F))
  W_vol   = 0.5 * K * (0.5 * (J^2 - 1) - log(J))
  W_dev   = 0.5 * G * (I_1_bar - 3.)
  return W_vol + W_dev
end

function energy(::NeoHookean, props::V, C::T) where {V <: AbstractArray{<:Number, 1}, T <: SymmetricTensor{2, 3, <:Number}}
  K, G    = props[1], props[2]
  I_3     = det(C)
  I_1_bar = tr(I_3^(-1. / 3.) * C)
  W_vol   = 0.5 * K * (0.5 * (I_3 - 1) - 0.5 * log(I_3))
  W_dev   = 0.5 * G * (I_1_bar - 3.)
  return W_vol + W_dev
end

function pk1_stress(::NeoHookean, props::V, F::T) where {V <: AbstractArray{<:Number, 1}, T <: Tensor{2, 3, <:Number, 9}}
  K, G    = props[1], props[2]
  J       = det(F)
  J_23    = J^(-2. / 3.)
  I_1     = tr(tdot(F))
  F_inv_T = inv(F)'
  P       = 0.5 * K * (J^2 - 1.) * F_inv_T + 
            G * J_23 * (F - (1. / 3.) * I_1 * F_inv_T)
  return P
end

function pk2_stress(::NeoHookean, props::V, C::T) where {V <: AbstractArray{<:Number, 1}, T <: SymmetricTensor{2, 3, <:Number, 6}}
  K, G     = props[1], props[2]
  I_1, I_3 = tr(C), det(C)
  J_23     = I_3^(-1. / 3.)
  C_inv    = inv(C)
  I        = one(SymmetricTensor{2, 3, eltype(C), 6})
  S        = 0.5 * K * (I_3 - 1.) * C_inv + 
             G * J_23 * (I - (1. / 3.) * I_1 * C_inv)
  return S
end

function cauchy_stress(::NeoHookean, props::V, F::T) where {V <: AbstractArray{<:Number, 1}, T <: Tensor{2, 3, <:Number, 9}}
  K, G = props[1], props[2]
  J    = det(F)
  J_53 = J^(-5. / 3.)
  B    = dott(F)
  I_1  = tr(B)
  I    = one(SymmetricTensor{2, 3, eltype(F), 6})
  σ    = 0.5 * K * (J - 1. / J) * I + 
         G * J_53 * (B - (1. / 3.) * I_1 * I)
  return σ
end