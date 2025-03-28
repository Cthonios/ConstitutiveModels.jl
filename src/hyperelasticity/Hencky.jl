struct Hencky <: AbstractHyperelasticity{2}
end

property_map(::Hencky) = Symbol[
  Symbol("bulk modulus"),
  Symbol("shear modulus")
]

@inline function helmholtz_free_energy(::Hencky, props, ∇u)
  K, G  = props[1], props[2]
  F     = deformation_gradient(∇u)
  B     = left_cauchy_green(F)
  E     = hencky_strain(B)
  # E     = 0.5 * log(tdot(F.val))
  # J     = det(F)
  trE   = tr(E)
  E_dev = E - (1 / 3) * trE * one(E)
  ψ_vol = 0.5 * K * trE^2
  ψ_dev = G * dcontract(E_dev, E_dev)
  return ψ_vol + ψ_dev
end

@inline function pk1_stress(::Hencky, props, ∇u)
  K, G  = props[1], props[2]
  F     = deformation_gradient(∇u)
  B     = left_cauchy_green(F)
  E     = hencky_strain(B)
  trE   = tr(E)
  E_dev = E - (1 / 3) * trE * one(E)
  σ = K * trE * one(E) + 2 * G * E_dev
  J = det(F)
  P = J * dot(σ, inv(F)')
  return P
end

# @inline function cauchy_stress(::Hencky, props, E)
#   K, G  = props[1], props[2]
#   trE   = tr(E)
#   E_dev = E - (1 / 3) * trE * one(E)
#   σ = K * trE * one(E) + 2 * G * E_dev
#   return σ
# end