@muladd begin

struct NeoHookean <: AbstractHyperelasticity{2}
end

property_map(::NeoHookean) = Symbol[
  Symbol("bulk modulus"),
  Symbol("shear modulus")
]

@inline function helmholtz_free_energy(::NeoHookean, props, ∇u)
  K, G    = props[1], props[2]
  F       = deformation_gradient(∇u)
  C       = right_cauchy_green(F)
  J       = det(F)
  I_1_bar = tr(J^(-2. / 3.) * C)
  W_vol   = 0.5 * K * (0.5 * (J^2 - 1) - log(J))
  W_dev   = 0.5 * G * (I_1_bar - 3.)
  ψ       = W_vol + W_dev
  return ψ
end

@inline function pk1_stress(::NeoHookean, props, ∇u)
  K, G    = props[1], props[2]
  F       = deformation_gradient(∇u)
  C       = right_cauchy_green(F)
  J       = det(F)
  J_23    = J^(-2. / 3.)
  I_1     = tr(C)
  F_inv_T = inv(F)'
  P       = 0.5 * K * (J^2 - 1.) * F_inv_T + 
            G * J_23 * (F - (1. / 3.) * I_1 * F_inv_T)
  return P
end

# TODO make efficient cauchy stress/pk2 stress implementation

end