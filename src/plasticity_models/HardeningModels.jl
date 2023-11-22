# """
# """
# abstract type HardeningModel{NProps, NStateVars} <: PlasticityModel{NProps, NStateVars} end
"""
"""
abstract type IsotropicHardeningModel{NProps, NStateVars} <: HardeningModel{NProps, NStateVars} end
function radius end
function slope end

function objective(
  model::Mod, 
  props, μ, σ_eff, α_old,
  u, p
) where Mod <: IsotropicHardeningModel

  Δγ = u[1]
  α_new = α_old + sqrt(2. / 3.) * Δγ
  g = σ_eff - radius(model, props, α_new) - 2. * μ * Δγ
  return [g]
end

function update(
  model::Mod, yield_surf::J2YieldSurface, 
  props, μ, σ, α_old
) where Mod <: IsotropicHardeningModel
  
  σ_eff = effective_stress(yield_surf, σ)
  
  x0 = SVector{1, eltype(σ_eff)}(0.0)
  f = (u, p) -> objective(model, props, μ, σ_eff, α_old, u, p)
  problem = NonlinearProblem(f, x0)
  sol = solve(problem, NewtonRaphson())
  Δγ = sol.u[1]
  return Δγ
end

"""
"""
struct NoIsotropicHardening <: IsotropicHardeningModel{0, 0}
end

"""
"""
function NoIsotropicHardening(::D) where D <: Dict
  return NoIsotropicHardening(), SVector{0, Float64}(), SVector{0, Float64}()
end

"""
"""
radius(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray = sqrt(2. / 3.) * props[1]
"""
"""
slope(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray  = 0.0
"""
Specialized method for update so no non-linear solve takes place
"""
function update(
  ::NoIsotropicHardening, yield_surf::J2YieldSurface, 
  props, μ, σ, α_old
)

  f = yield_surface(yield_surf, props, σ, α_old)

  if f <= 0.0
    Δγ = 0.0
  else
    Δγ = f / (2. * μ)
  end
  return Δγ
end


"""
"""
struct LinearIsotropicHardening <: IsotropicHardeningModel{1, 0}
end

"""
"""
function LinearIsotropicHardening(inputs::D) where D <: Dict
  @assert "hardening modulus" in keys(inputs)
  H = inputs["hardening modulus"]
  return LinearIsotropicHardening(), SVector{1, Float64}(H), SVector{0, Float64}()
end

"""
"""
radius(::LinearIsotropicHardening, props, eqps) = sqrt(2. / 3.) * (props[1] + props[2] * eqps)
"""
"""
slope(::LinearIsotropicHardening, props, eqps)  = props[2]

"""
"""
function update(
  ::LinearIsotropicHardening, yield_surf::J2YieldSurface, 
  props, μ, σ, α_old
)
  H = props[2]
  f = yield_surface(yield_surf, props, σ, α_old)

  if f <= 0.0
    Δγ = 0.0
  else
    Δγ = f / (2. * μ * (1. + (H / (3. * μ))))
  end
  return Δγ
end

struct VoceIsotropicHardening <: IsotropicHardeningModel{2, 0}
end

function VoceIsotropicHardening(inputs::D) where D <: Dict
  @assert "A" in keys(inputs)
  @assert "n" in keys(inputs)
  A = inputs["A"]
  n = inputs["n"]

  return VoceIsotropicHardening(), SVector{2, Float64}((A, n)), SVector{0, Float64}()
end

radius(::VoceIsotropicHardening, props, eqps) = sqrt(2. / 3.) * (props[1] + props[2] * (1. - exp(-props[3] * eqps)))
slope(::VoceIsotropicHardening, props, eqps)  = props[2] * props[3] * exp(-props[3] * eqps)