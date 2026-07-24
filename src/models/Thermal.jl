struct Thermal{C, F} <: AbstractHyperelasticModel
    heat_capacity_model::C
    heat_flux_model::F
end

function initialize_props(model::Thermal, inputs::Dict{String})
    return [
        inputs["density"],
        # inputs["reference temperature"],
        # inputs["specific heat capacity"],
        initialize_props(model.heat_capacity_model, inputs)...,
        initialize_props(model.heat_flux_model, inputs)...
    ]
end

num_properties(model::Thermal) = num_properties(model.heat_capacity_model) + num_properties(model.heat_flux_model) + 1
num_state_variables(model::Thermal) = 0

function dissipation(
    model::Thermal,
    props, Z_old, Z_new, Δt, ∇u, θ, ∇θ
)
    props = module_props(model.heat_flux_model, props, 4)
    return dissipation(model.heat_flux_model, props, ∇u, θ, ∇θ)
end

function entropy(
    model::Thermal,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    # θ_0, c = props[2], props[3]
    # η = (c / θ_0) * (θ - θ_0)
    # return η
    props = module_props(model.heat_capacity_model, props, 2)
    return entropy(model.heat_capacity_model, props, ∇u, θ)
end

function heat_capacity(
    model::Thermal,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    # θ_0, c = props[2], props[3]
    # return c * (θ / θ_0)
    props = module_props(model.heat_capacity, props, 2)
    return heat_capacity(model.heat_capacity_model, props, ∇u, θ)
end

function heat_flux(
    model::Thermal,
    props, Z_old, Z_new, Δt, ∇u, θ, ∇θ
)
    props = module_props(model.heat_flux_model, props, 4)
    return heat_flux(model.heat_flux_model, props, ∇u, θ, ∇θ)
end

function helmholtz_free_energy(
    model::Thermal,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    # θ_0, c = props[2], props[3]
    # ψ = -(c / (2 * θ_0)) * (θ - θ_0)^2
    # return ψ
    props = module_props(model.heat_capacity_model, props, 2)
    return helmholtz_free_energy(model.heat_capacity_model, props, ∇u, θ)
end
