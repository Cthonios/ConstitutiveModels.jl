abstract type AbstractHyperelasticity{NP} <: AbstractMechanicalModel{NP, 0} end

function initialize_properties(model::AbstractHyperelasticity{NP}, inputs) where NP
  prop_map = property_map(model)
  props = SVector{NP, Float64}([inputs[prop_map[n]] for n in axes(prop_map, 1)])
  return props
end

function cauchy_stress(model::AbstractHyperelasticity, props, ε::SymmetricTensor)
  # F = deformation_gradient(∇u)
  # J = det(F)
  # P = pk1_stress(model, props, ∇u)
  # return (1 / J) * P * F'
  return Tensors.gradient(z -> helmholtz_free_energy(model, props, z), ε)
end

function cauchy_stress(model::AbstractHyperelasticity, props, ∇u)
  F = ∇u + one(∇u)
  P = pk1_stress(model, props, ∇u)
  J = det(F)
  return (1 / J) * dot(P, F')
end

function pk1_stress(model::AbstractHyperelasticity, props, ∇u)
  return Tensors.gradient(z -> helmholtz_free_energy(model, props, z), ∇u)
end

include("Hencky.jl")
include("LinearElastic.jl")
include("NeoHookean.jl")
