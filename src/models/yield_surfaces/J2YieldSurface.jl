struct J2YieldSurface <: YieldSurface{1, 7}
end

function J2YieldSurface(inputs::D, type::Type{ArrType}) where {D <: Dict, ArrType}
  yield_stress = inputs["yield stress"]
  model = J2YieldSurface()
  props = initialize_props((yield_stress,), type)
  state = initialize_state(model, type)
  return model, props, state
end

function effective_stress(::J2YieldSurface, σ::M) where M <: AbstractArray{<:Number, 2}
  σ_vm = norm(dev(σ))
  return σ_vm
end

# function yield_surface(model::J2YieldSurface, props::V1, σ::V2, α_old) where {V1, V2}
#   σ_eff = effective_stress(model, σ)
#   f = σ_eff # TODO add hardening
# end