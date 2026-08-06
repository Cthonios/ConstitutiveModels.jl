struct Hyperelastic{
    M <: AbstractHyperelasticity
} <: AbstractHyperelasticModel
    model::M
end

function initialize_props(model::Hyperelastic, inputs::Dict{String})
    return [get_property(inputs, "density"), initialize_props(model.model, inputs)...]
end

num_properties(model::Hyperelastic) = 1 + num_properties(model.model)
num_state_variables(model::Hyperelastic) = 0

function property_names(model::Hyperelastic)
    return ["density", property_names(model.model)...]
end

function cauchy_stress(
    model::Hyperelastic,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    props = module_props(model.model, props, 2)
    return cauchy_stress(model.model, props, ∇u, θ)
end

function helmholtz_free_energy(
    model::Hyperelastic,
    props, Z_old, Z_new, Δt, ∇u, θ
)   props = module_props(model.model, props, 2)
    return helmholtz_free_energy(model.model, props, ∇u, θ)
end

function pk1_stress(
    model::Hyperelastic,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    props = module_props(model.model, props, 2)
    return pk1_stress(model.model, props, ∇u, θ)
end

function material_tangent(
    model::Hyperelastic,
    props, Z_old, Z_new, Δt, ∇u, θ
)    
    props = module_props(model.model, props, 2)
    return material_tangent(model.model, props, ∇u, θ)
end

function p_wave_modulus(model::Hyperelastic, props)
    props = module_props(model.model, props, 2)
    return p_wave_modulus(model.model, props)
end
