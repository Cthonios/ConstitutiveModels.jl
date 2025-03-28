@muladd begin

struct LinearElastic <: AbstractHyperelasticity{2}
end

property_map(::LinearElastic) = Symbol[
  Symbol("youngs modulus"),
  Symbol("poissons ratio")
]

function initialize_properties(model::LinearElastic, inputs)
  prop_map = property_map(model)
  E = inputs[prop_map[1]]
  ν = inputs[prop_map[2]]
  μ = E / (2 * (1 + ν))
  λ = E * ν / ((1 + ν) * (1 - 2 * ν))
  return SVector{2, Float64}(λ, μ)
end

function helmholtz_free_energy(model::LinearElastic, props, ∇u::Tensor)
  ε = linear_strain(∇u)
  return helmholtz_free_energy(model, props, ε)
end

function helmholtz_free_energy(model::LinearElastic, props, ε::SymmetricTensor)
  # λ, μ = props[1], props[2]
  # ψ    = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)
  ψ = deviatoric_helmholtz_free_energy(model, props, ε) + 
      volumetric_helmholtz_free_energy(model, props, ε)
  return ψ
end

function deviatoric_helmholtz_free_energy(::LinearElastic, props, ε::SymmetricTensor)
  μ = props[2]
  ψ = μ * dcontract(ε, ε)
  return ψ
end

function volumetric_helmholtz_free_energy(::LinearElastic, props, ε::SymmetricTensor)
  λ = props[1]
  ψ = 0.5 * λ * tr(ε)^2
  return ψ
end

end