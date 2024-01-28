struct VoceIsotropicHardening <: IsotropicHardening{2, 0}
end

function VoceIsotropicHardening(inputs::D, type::Type{ArrType}) where {D <: Dict, ArrType}
  A = inputs["A"]
  n = inputs["n"]
  model = VoceIsotropicHardening()
  props = initialize_props((A, n), type)
  state = initialize_state(model, type)
  return model, props, state
end

energy(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray = props[2] * eqps + (props[2] - props[1]) * props[3] * exp(-eqps / props[3])
radius(::VoceIsotropicHardening, props, eqps) = sqrt(2. / 3.) * (props[1] + props[2] * (1. - exp(-props[3] * eqps)))
slope(::VoceIsotropicHardening, props, eqps)  = sqrt(2. / 3.) * props[2] * props[3] * exp(-props[3] * eqps)
