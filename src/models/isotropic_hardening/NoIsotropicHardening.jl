struct NoIsotropicHardening <: IsotropicHardening{0, 0}
end

function NoIsotropicHardening(::D, type::Type{ArrType}) where {D <: Dict, ArrType}
  model = NoIsotropicHardening()
  props = initialize_props((), type)
  state = initialize_state(model, type)
  return model, props, state
end

energy(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray = props[1] * eqps
radius(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray = sqrt(2. / 3.) * props[1]
slope(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray  = 0.0

function hardening_increment(
  model::NoIsotropicHardening,
  props, μ, σ_eff::T, α_old
) where T <: Number

  f = σ_eff - radius(model, props, α_old)

  if f <= zero(T)
    Δγ = 0.0
  else
    Δγ = f / (2. * μ)
  end
  return Δγ
end
