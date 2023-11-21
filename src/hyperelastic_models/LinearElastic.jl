struct LinearElastic <: HyperelasticModel{2, 0}
end

function LinearElastic(inputs::D) where D <: Dict
  @assert "youngs modulus" in keys(inputs)
  @assert "poissons ratio" in keys(inputs)

  E = inputs["youngs modulus"]
  ν = inputs["poissons ratio"]
  
  λ = E * ν / ((1. + ν) * (1. - 2. * ν)) 
  μ = E / (2. * (1. + ν))

  # TODO add props check
  props = @SVector [
    λ,
    μ
  ]
  return LinearElastic(), props
end

function energy(model::LinearElastic, props::V, F::Tensor{2, 3, T, 9}) where {V <: AbstractArray, T <: Number}
  # λ, μ = props[1], props[2]
  # I = one(SymmetricTensor{2, 3, eltype(F), 6})
  I = one(Tensor{2, 3, eltype(F), 9})
  ∇u = F - I
  # ε = symmetric(∇u)
  ε = 0.5 * (∇u + ∇u')
  # trε  = tr(ε)
  # # return energy(model, props, ε)
  # W = 0.5 * λ * trε^2 + μ * dcontract(ε, ε)
  # return W
  return energy_enzyme_specialization(model, props, ε)
end

# shoudl be general AD energy method
# function energy(::LinearElastic, props::V, ε::SymmetricTensor{2, 3, T, 6}) where {V <: AbstractArray, T <: Number}
#   λ, μ = props[1], props[2]
#   trε  = tr(ε)
#   W = 0.5 * λ * trε^2 + μ * dcontract(ε, ε)
#   return W
# end

# Specialization for enzyme since enzyme double counts some things in
# symmetric gradients. May need a custom rule 
function energy_enzyme_specialization(::LinearElastic, props::V, ε::Tensor{2, 3, T, 9}) where {V <: AbstractArray, T <: Number}
  λ, μ = props[1], props[2]
  trε  = tr(ε)
  W = 0.5 * λ * trε^2 + μ * dcontract(ε, ε)
  return W
end

function cauchy_stress(model::LinearElastic, props::V, F::T) where {V <: AbstractArray, T <: Tensor{2, 3, <:Number, 9}}
  I = one(SymmetricTensor{2, 3, eltype(F), 6})
  ∇u = F - I
  ε = symmetric(∇u)
  # ε = 0.5 * (∇u + ∇u')
  # return cauchy_stress(model, props, ε)
  λ, μ = props[1], props[2]
  I = one(SymmetricTensor{2, 3, eltype(ε), 6})
  trε  = tr(ε)
  σ = λ * trε * I + 2. * μ * ε
end

function cauchy_stress(::LinearElastic, props::V, ε::T) where {V <: AbstractArray, T <: SymmetricTensor{2, 3, <:Number, 6}}
  λ, μ = props[1], props[2]
  I = one(SymmetricTensor{2, 3, eltype(ε), 6})
  trε  = tr(ε)
  σ = λ * trε * I + 2. * μ * ε
  return σ
end

# making it the same as cauchy stress for the linear case
function pk1_stress(model::LinearElastic, props::V, F::T) where {V <: AbstractArray, T <: Tensor{2, 3, <:Number, 9}}
  return cauchy_stress(model, props, F)
end
