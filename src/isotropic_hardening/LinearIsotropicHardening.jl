struct LinearIsotropicHardening <: IsotropicHardening{1, 0}
end

function LinearIsotropicHardening(_)
  return LinearIsotropicHardening()
end

function initialize_properties(::LinearIsotropicHardening, inputs::D) where D <: Dict{Symbol, Any}
  H = inputs[Symbol("hardening modulus")]
  return SVector{1, Float64}(H)
end

energy(::LinearIsotropicHardening, props::V, eqps) where V <: AbstractArray = props[1] * eqps + 0.5 * props[2] * eqps * eqps
radius(::LinearIsotropicHardening, props, eqps) = sqrt(2. / 3.) * (props[1] + props[2] * eqps)
slope(::LinearIsotropicHardening, props, eqps)  = props[2]

function hardening_increment(
  model::LinearIsotropicHardening,
  props, μ, σ_eff::T, α_old
) where T <: Number

  H = props[2]

  f = σ_eff - radius(model, props, α_old)

  if f <= zero(T)
    Δγ = 0.0
  else
    Δγ = f / (2. * μ * (1. + (H / (3. * μ))))
  end
  return Δγ
end
