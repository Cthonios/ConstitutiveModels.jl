abstract type ThermalModel{NP, NS} <: ConstitutiveModel{NP, NS} end

function heat_flux end

struct FouriersLaw <: ThermalModel{1, 0}
end

function initialize_props(::FouriersLaw, inputs::D) where {D <: Dict{Symbol, Any}}
  k = inputs[Symbol("thermal conductivity")]
  return initialize_props((k,))
end

function heat_flux(::FouriersLaw, props, Δt, T, ∇T, state_old)

  # unpack props
  k = props[1]

  # caclulate heat flux
  q = -k * ∇T

  # dummy state
  state_new = typeof(state_old)()

  return q, state_new
end
