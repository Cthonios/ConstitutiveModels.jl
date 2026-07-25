struct LinearThermoelastic{
    E <: LinearElasticity,
    C <: SimpleHeatCapacity,
    H <: FouriersLaw,
    T <: LinearThermalExpansion
} <: AbstractLinearElasticModel
    elastic_model::E
    heat_capacity_model::C
    heat_flux_model::H
    thermal_expansion_model::T
end

function LinearThermoelastic(
    ;
    elasticity_symmetry = nothing,
    heat_conduction_symmetry = Isotropy{2}(),
    thermal_expansion_symmetry = Isotropy{2}()
)
    return LinearThermoelastic(
        LinearElasticity(),
        SimpleHeatCapacity(),
        FouriersLaw(heat_conduction_symmetry),
        LinearThermalExpansion(thermal_expansion_symmetry)
    )
end

function initialize_props(::LinearThermoelastic, inputs::Dict{String})
    # elastic_props = ElasticConstants(inputs)
    # α = inputs["coefficient of thermal expansion"]
    # θ_0 = inputs["reference temperature"]
    # c = inputs["specific heat capacity"]
    # k = inputs["thermal conductivity"]
    # β = -3. * elastic_props.κ * α
    # return [
    #     inputs["density"],
    #     elastic_props.λ, elastic_props.μ,
    #     β, θ_0, c, k
    # ]
    return [
        get_property(inputs, "density"),
        initialize_props(model.heat_capacity_model, inputs)...,
        initialize_props(model.heat_flux_model, inputs)...,
        initialize_props(model.thermal_expansion, inputs)...,
        initialize_props(model.elastic_model, inputs)...,
    ]
end

num_properties(::LinearThermoelastic) = 1 + num_properties(model.elastic_model) +
    num_properties(model.heat_capacity_model) + num_properties(model.heat_flux_model) +
    num_prony_terms(model.thermal_expansion_model)
num_state_variables(::LinearThermoelastic) = 0

function dissipation(
    ::LinearThermoelastic,
    props, Z_old, Z_new, Δt, ε, θ, ∇θ
)
    props = module_props(model.heat_flux_model, props, 7)
    return dissipation(model.heat_flux_model, props, ε, θ, ∇θ)
end

function heat_capacity(
    model::LinearThermoelastic,
    props, Z_old, Z_new, Δt, ε, θ, ∇θ
)

end

function heat_flux(
    model::LinearThermoelastic,
    props, Z_old, Z_new, Δt, ε, θ, ∇θ
)
    props_k = module_props(model.heat_flux_model, props, 7)
    return heat_flux(model.heat_flux_model, props_k, ε, θ, ∇θ)
end

function helmholtz_free_energy(
    model::LinearThermoelastic,
    props, Z_old, Z_new, Δt, ε, θ
)
    _, λ, μ, β, θ_0, c, _ = props
    props_m = SVector{2, typeof(λ)}(λ, μ)
    ψ_m = helmholtz_free_energy(model.elastic_model, props_m, ε, θ)
    ψ_t = -(c / (2. * θ_0)) * (θ - θ_0)^2
    ψ_mixed = β * (θ - θ_0) * tr(ε)
    return ψ_m + ψ_t + ψ_mixed
end

function p_wave_modulus(model::LinearThermoelastic, props)
    props = module_props(model.elastic_model, props, 2)
    return p_wave_modulus(model.elastic_model, props)
end