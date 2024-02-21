abstract type IsotropicHardening{NP, NS} <: ConstitutiveModel{NP, NS} end
function radius end
function slope end

function hardening_objective(
  model::Mod,
  props, μ, σ_eff, α_old, 
  Δγ
) where Mod <: IsotropicHardening

  α_new = α_old + sqrt(2. / 3.) * Δγ
  g = σ_eff - radius(model, props, α_new) - 2. * μ * Δγ
  return g
end

function hardening_slope(
  model::Mod,
  props, μ, σ_eff, α_old, 
  Δγ
) where Mod <: IsotropicHardening
  α_new = α_old + sqrt(2. / 3.) * Δγ
  return -slope(model, props, Δγ) - 2. * μ
  # return -2. * μ * Δγ
end

function hardening_increment(
  model::IsotropicHardening,
  props, μ, σ_eff::T, α_old
) where T <: Number
  
  f = x -> hardening_objective(model, props, μ, σ_eff, α_old, x)
  g = x -> hardening_slope(model, props, μ, σ_eff, α_old, x)
  Δγ = solve_hardening(f, g, α_old)
  return Δγ
end

include("LinearIsotropicHardening.jl")
include("NoIsotropicHardening.jl")
include("SwiftVoceIsotropicHardening.jl")
include("VoceIsotropicHardening.jl")
