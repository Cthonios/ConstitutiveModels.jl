struct Hencky <: HyperelasticModel{2, 0}
end

function Hencky(_)
  return Hencky()
end

function initialize_props(::Hencky, inputs::D) where {D <: Dict{Symbol, Any}}
  K = inputs[Symbol("bulk modulus")]
  G = inputs[Symbol("shear modulus")]
  # return initialize_props((K, G), type)
  return SVector{2, Float64}((K, G))
end

function helmholtz_free_energy(
  ::Hencky,
  props::V1, Δt, F::M, θ, state_old
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2},        
}
  K, G = props[1], props[2]
  # I = one(typeof(F))
  J = det(F)
  trE_e = log(J)
  E_dev = 0.5 * J^(-2. / 3.) * log(tdot(F))
  # E     = E_dev + (1. / 3.) * trE_e * I

  ψ_vol = 0.5 * K * trE_e^2
  ψ_dev = G * dcontract(E_dev, E_dev)

  # dummy state
  state_new = SVector{0, eltype(state_old)}()

  return ψ_vol + ψ_dev, state_new
end

function helmholtz_free_energy_old(
  ::Hencky,
  props::V1, Δt, F::M, θ, state_old
) where {
  V1 <: AbstractArray{<:Number, 1},
  M  <: AbstractArray{<:Number, 2},        
}
  K, G = props[1], props[2]
  I = one(typeof(F))
  J = det(F)
  C_bar = tdot(J^(-2. / 3.) * F)
  E_bar = 0.5 * log(C_bar)

  # ψ = λ * tr(E)^2 + 2. * μ * dcontract(E, E)
  # ψ = K * tr(E) * I + 2. * G * dev(E)
  ψ = 0.5 * K * tr(E_bar)^2 + G * dcontract(dev(E), dev(E))

  # dummy state
  state_new = SVector{0, Float64}()

  return ψ, state_new
  # ∇u = F - I
  # ε = 0.5 * (∇u + ∇u')
  # W = 0.5 * λ * tr(ε)^2 + μ * dcontract(ε, ε)
  # return W
end