struct NeoHookean <: HyperelasticModel{2, 0}
end

function read_properties(inputs::D) where D <: Dict{String, Any}
  bulk_modulus  = inputs["bulk modulus"]
  shear_modulus = inputs["shear modulus"]
  return bulk_modulus, shear_modulus
end

function read_properties(inputs::D) where D <: Dict{Symbol, Any}
  bulk_modulus  = inputs[Symbol("bulk modulus")]
  shear_modulus = inputs[Symbol("shear modulus")]
  return bulk_modulus, shear_modulus
end

function NeoHookean(inputs::D) where D <: Dict
  # @assert "bulk modulus" in keys(inputs)
  # @assert "shear modulus" in keys(inputs)
  # bulk_modulus  = inputs["bulk modulus"]
  # shear_modulus = inputs["shear modulus"]
  bulk_modulus, shear_modulus = read_properties(inputs)

  model = NeoHookean()
  props = initialize_properties(model, [
    bulk_modulus, 
    shear_modulus
  ])
    # TODO add prop checks
  state = initialize_state(model) # make type parametric TODO
  
  return model, props, state
end
 
# pure method for energy
function energy(
  ::NeoHookean, 
  props::V1, 
  F::M, 
  state::V2
) where {V1 <: AbstractArray{<:Number, 1}, 
         M  <: AbstractArray{<:Number, 2}, 
         V2 <: AbstractArray{<:Number, 1}}

  K, G    = props[1], props[2]
  J       = det(F)
  I_1_bar = tr(J^(-2. / 3.) * tdot(F))
  W_vol   = 0.5 * K * (0.5 * (J^2 - 1) - log(J))
  W_dev   = 0.5 * G * (I_1_bar - 3.)
  return W_vol + W_dev
end

function pk1_stress(
  ::NeoHookean, 
  props::V1, 
  F::M, 
  state::V2
) where {V1 <: AbstractArray{<:Number, 1}, 
         M  <: AbstractArray{<:Number, 2},
         V2 <: AbstractArray{<:Number, 1}}

  K, G    = props[1], props[2]
  J       = det(F)
  J_23    = J^(-2. / 3.)
  I_1     = tr(tdot(F))
  F_inv_T = inv(F)'
  P       = 0.5 * K * (J^2 - 1.) * F_inv_T + 
            G * J_23 * (F - (1. / 3.) * I_1 * F_inv_T)
  return P
end

function cauchy_stress(
  ::NeoHookean, 
  props::V1, 
  F::M, 
  state::V2
) where {V1 <: AbstractArray{<:Number, 1}, 
         M  <: AbstractArray{<:Number, 2},
         V2 <: AbstractArray{<:Number, 1}}

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
