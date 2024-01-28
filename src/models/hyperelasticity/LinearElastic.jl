struct LinearElastic <: HyperelasticModel{2, 0}
end

function LinearElastic(inputs::D, type::Type{ArrType}) where {D <: Dict, ArrType}
  E = inputs["youngs modulus"]
  ν = inputs["poissons ratio"]
  
  λ = E * ν / ((1. + ν) * (1. - 2. * ν)) 
  μ = E / (2. * (1. + ν))

  # TODO add props check
  model = LinearElastic()
  props = initialize_props((λ, μ), type)
  state = initialize_state(model, type)
  return model, props, state
end

function helmholtz_free_energy(
  ::LinearElastic,
  props::V1, F::M
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2},        
}
  λ, μ = props[1], props[2]
  I = one(typeof(F))
  ∇u = F - I
  ε = 0.5 * (∇u + ∇u')
  W = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)
  return W
end

function cauchy_stress(
  ::LinearElastic,
  props::V1, F::M
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2}
}
  λ, μ = props[1], props[2]
  I = one(typeof(F))
  ∇u = F - I
  ε = 0.5 * (∇u + ∇u')
  σ = λ * tr(ε) * I + 2. * μ * ε
  return σ
end

function cauchy_stress(
  ::LinearElastic,
  props::V1, ε::SymmetricTensor{2, 3, T, 6}
) where {
  V1 <: AbstractArray{<:Number, 1},
  T  <: Number
}
  λ, μ = props[1], props[2]
  I = one(SymmetricTensor{2, 3, T, 6})
  σ = λ * tr(ε) * I + 2. * μ * ε
  return σ
end

function pk1_stress(
  model::LinearElastic,
  props::V1, F::M
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2}
}
  return cauchy_stress(model, props, F)
end
