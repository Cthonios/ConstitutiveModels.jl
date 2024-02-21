struct SwiftVoceIsotropicHardening <: IsotropicHardening{3, 0}
end

function SwiftVoceIsotropicHardening(inputs::D, type::Type{ArrType}) where {D <: Dict, ArrType}
  H = inputs["hardening modulus"]
  A = inputs["A"]
  n = inputs["n"]
  model = SwiftVoceIsotropicHardening()
  props = initialize_props((H, A, n), type)
  state = initialize_state(model, type)
  return model, props, state
end

radius(::SwiftVoceIsotropicHardening, props, eqps) = sqrt(2. / 3.) * (props[1] + props[2] * eqps + props[3] * (1. - exp(-props[4] * eqps)))
slope(::SwiftVoceIsotropicHardening, props, eqps)  = sqrt(2. / 3.) * (props[2] + props[3] * props[4] * exp(-props[4] * eqps))
