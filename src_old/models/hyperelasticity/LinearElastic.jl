struct LinearElastic <: HyperelasticModel{2, 0}
end

function LinearElastic(_)
  return LinearElastic()
end

function initialize_props(::LinearElastic, inputs::D) where {D <: Dict{Symbol, Any}}
  E = inputs[Symbol("youngs modulus")]
  ν = inputs[Symbol("poissons ratio")]
  λ = E * ν / ((1. + ν) * (1. - 2. * ν)) 
  μ = E / (2. * (1. + ν))
  return SVector{2, Float64}((λ, μ))
end

function helmholtz_free_energy(
  ::LinearElastic,
  props::V1, Δt, F::M, θ, state_old::V2
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2},  
  V2 <: AbstractArray{<:Number, 1}       
}
  # unpack properties
  λ, μ = props[1], props[2]

  # kinematics
  I = one(typeof(F))
  ∇u = F - I
  # ε = 0.5 * (∇u + ∇u')
  ε = symmetric(∇u)

  # constitutive
  ψ = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)
  # ψ = 0.5 * λ * tr(ε)^2

  # dummy state
  state_new = V2()

  return ψ, state_new
end

# function pk1_stress(
#   ::LinearElastic,
#   props::V1, Δt, F::M, θ, state_old::V2
# ) where {
#   V1 <: AbstractArray{<:Number, 1},
#   M  <: AbstractArray{<:Number, 2},  
#   V2 <: AbstractArray{<:Number, 1}       
# }
#   # unpack properties
#   λ, μ = props[1], props[2]

#   # kinematics
#   I = one(typeof(F))
#   ∇u = F - I
#   # ε = symmetric(∇u)
#   ε = 0.5 * (∇u + ∇u')

#   # constitutive
#   σ = λ * tr(ε) * I + 2. * μ * ε

#   # dummy state
#   state_new = V2()

#   return σ, state_new
# end

function cauchy_stress(
  ::LinearElastic,
  props::V1, Δt, F::M, θ, state_old::V2
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2},  
  V2 <: AbstractArray{<:Number, 1}       
}

  # unpack properties
  λ, μ = props[1], props[2]

  # kinematics
  I = one(typeof(F))
  ∇u = F - I
  ε = 0.5 * (∇u + ∇u')

  # constitutive
  σ = λ * tr(ε) * I + 2. * μ * ε

  # dummy state
  state_new = V2()

  return σ, state_new
end

# mainly used for plasticity models
function helmholtz_free_energy(
  ::LinearElastic,
  props::V1, ε::SymmetricTensor{2, 3, <:Number, 6}
) where {
  V1 <: AbstractArray{<:Number, 1}     
}
  # unpack properties
  λ, μ = props[1], props[2]

  # kinematics
  # I = one(typeof(ε))

  # constitutive
  ψ = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)

  return ψ
end

# mainly used for plasticity models
function cauchy_stress(
  ::LinearElastic,
  props, ε
)

  # unpack properties
  λ, μ = props[1], props[2]

  # constitutive
  I = one(SymmetricTensor{2, 3, eltype(ε), 6})
  σ = λ * tr(ε) * I + 2. * μ * ε

  return σ
end