struct NeoHookean <: HyperelasticModel{2, 0}
end

function initialize_props(::NeoHookean, inputs::D, type::Type{ArrType}) where {ArrType, D <: Dict{String, Any}}
  K = inputs["bulk modulus"]
  G = inputs["shear modulus"]
  return initialize_props((K, G), type)
end

function initialize_props(::NeoHookean, inputs::D, type::Type{ArrType}) where {ArrType, D <: Dict{Symbol, Any}}
  K = inputs[Symbol("bulk modulus")]
  G = inputs[Symbol("shear modulus")]
  return initialize_props((K, G), type)
end

# function helmholtz_free_energy(
function helmholtz_free_energy(
  ::NeoHookean,
  props::V1, Δt, F::M, θ, state_old::V2
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2},
  V2 <: AbstractArray{<:Number, 1}        
}

  K, G    = props[1], props[2]
  J       = det(F)
  I_1_bar = tr(J^(-2. / 3.) * tdot(F))
  W_vol   = 0.5 * K * (0.5 * (J^2 - 1) - log(J))
  W_dev   = 0.5 * G * (I_1_bar - 3.)
  ψ       = W_vol + W_dev

  # dummy state
  state_new = V2()

  return ψ, state_new
end

function pk1_stress(
  ::NeoHookean, 
  props::V1, Δt, F::M, θ, state_old::V2
) where {
  V1 <: AbstractArray{<:Number, 1}, 
  M  <: AbstractArray{<:Number, 2},
  V2 <: AbstractArray{<:Number, 1}
}

  K, G    = props[1], props[2]
  J       = det(F)
  J_23    = J^(-2. / 3.)
  I_1     = tr(tdot(F))
  F_inv_T = inv(F)'
  P       = 0.5 * K * (J^2 - 1.) * F_inv_T + 
            G * J_23 * (F - (1. / 3.) * I_1 * F_inv_T)

  # dummy state
  state_new = V2()

  return P, state_new
end

# TODO add analytic tangent
