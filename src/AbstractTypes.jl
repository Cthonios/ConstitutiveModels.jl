"""
$(TYPEDEF)
Base type in the package.
"""
abstract type ConstitutiveModel{NP, NS} end

function ConstitutiveModel(model::Type{<:ConstitutiveModel}, inputs::Dict{String, Any})
  new_inputs = Dict{Symbol, Any}()
  for (key, val) in pairs(inputs)
    new_inputs[Symbol(key)] = val
  end
  return ConstitutiveModel(model, new_inputs)
end

function ConstitutiveModel(model::Type{<:ConstitutiveModel}, inputs::Dict{Symbol, Any})
  return model(inputs)
end

num_properties(::ConstitutiveModel{NP, NS}) where {NP, NS} = NP
num_state_vars(::ConstitutiveModel{NP, NS}) where {NP, NS} = NS

function initialize_properties(model::T, inputs::Dict{String, Any}) where T <: ConstitutiveModel
  new_inputs = Dict{Symbol, Any}()
  for (key, val) in pairs(inputs)
    new_inputs[Symbol(key)] = val
  end
  # return new_inputs
  return initialize_properties(model, new_inputs)
end

"""
$(TYPEDEF)
Base type for purely mechanical models.
"""
abstract type AbstractMechanicalModel{NP, NS} <: ConstitutiveModel{NP, NS} end

function cauchy_stress(model::AbstractMechanicalModel, props, ∇u, θ, state, Δt)
  F        = deformation_gradient(∇u)
  J        = det(F)
  P, state = pk1_stress(model, props, ∇u, θ, state, Δt)
  σ        = (1 / J) * dot(P, F')
  return σ, state
end

"""
$(TYPEDEF)
Base type for purely thermal models.
"""
abstract type AbstractThermalModel{NP, NS} <: ConstitutiveModel{NP, NS} end
