struct Hyperelastic{NP, S} <: AbstractMechanicalModel{NP, 0}
  strain_energy::S
end

function Hyperelastic(strain_energy::AbstractHyperelasticity)
  NP = num_properties(strain_energy)
  return Hyperelastic{NP + 1, typeof(strain_energy)}(strain_energy)
end

function Hyperelastic(inputs::Dict{Symbol, Any})
  strain_energy = eval(Symbol(inputs[Symbol("strain energy density")]))()
  return Hyperelastic{num_properties(strain_energy) + 1, typeof(strain_energy)}(strain_energy)
end

property_map(::Hyperelastic) = Symbol[
  Symbol("density"),
  Symbol("strain energy density")
]

function initialize_properties(::Hyperelastic, inputs::Dict{Symbol, Any})
  strain_energy_input = Symbol(inputs[Symbol("strain energy density")])
  strain_energy = eval(strain_energy_input)()
  strain_energy_props = initialize_properties(strain_energy, inputs)
  return vcat(inputs[:density], strain_energy_props...)
end

initialize_state(::Hyperelastic) = SVector{0, Float64}()

@inline function helmholtz_free_energy(model::Hyperelastic, props, ∇u, θ, state, Δt)
  return @views helmholtz_free_energy(model.strain_energy, props, ∇u), state
end

@inline function pk1_stress(model::Hyperelastic, props, ∇u, θ, state, Δt)
  return @views pk1_stress(model.strain_energy, props[2:end], ∇u), state
end

# NEED PK1 method and tangent or it will always fall back to AD!