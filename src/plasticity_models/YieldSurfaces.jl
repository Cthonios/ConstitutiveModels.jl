abstract type YieldSurface{NProps, NStateVars} <: PlasticityModel{NProps, NStateVars} end

"""
Always assumes the yield stress is the first prop
"""
yield_stress(::Mod, props::V) where {Mod <: YieldSurface, V <: AbstractArray} = props[1]

# 1 prop for yield stress
# 2 state vars
struct J2YieldSurface{NProps, Hard <: HardeningModel} <: YieldSurface{NProps, 7}
  hardening_model::Hard
end

# initialize_state(::J2YieldSurface) = @SVector [0.0, 0.0]

function J2YieldSurface(inputs::D) where D <: Dict
  @assert "yield stress" in keys(inputs)
  yield_stress = inputs["yield stress"]

  if "isotropic hardening model" in keys(inputs)
    error("Not implented yet")
  else
    hardening_model, hardening_props, hardening_state = NoIsotropicHardening(inputs)
  end

  model = J2YieldSurface{1 + length(hardening_props), typeof(hardening_model)}(hardening_model)
  props = initialize_properties(model, vcat([yield_stress], hardening_props))
  state = @SVector [0., 0., 0., 0., 0., 0., 0.]
  state = vcat(state, hardening_state)
  return model, props, state
end

function effective_stress(::J2YieldSurface, σ::M) where M <: AbstractArray{<:Number, 2}
  σ_vm = norm(dev(σ))
  return σ_vm
end

function yield_surface(model::J2YieldSurface, props::V1, σ::V2, α_old) where {V1 <: AbstractArray, V2 <: AbstractArray{<:Number, 2}}
  # σ_y = yield_stress(model, props)
  σ_eff = effective_stress(model, σ)
  surf_radius = radius(model.hardening_model, props, α_old)
  f = σ_eff - surf_radius
  return f
end

# TODO maybe a better interface?
function update(model::J2YieldSurface, props::V1, μ, σ::V2, α_old) where {V1 <: AbstractArray, V2 <: AbstractArray}
  f = yield_surface(model, props, σ, α_old)
  if f <= 0.0
    Δγ = 0.0
  else
    Δγ = f / (2. * μ)
  end
  return Δγ
end
