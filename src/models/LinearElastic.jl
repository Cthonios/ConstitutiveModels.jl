struct LinearElastic <: AbstractLinearElasticModel
    model::LinearElasticity
end

function LinearElastic()
    return LinearElastic(LinearElasticity())
end

function initialize_props(model::LinearElastic, inputs::Dict{String})
    return [get_property(inputs, "density"), initialize_props(model.model, inputs)...]
end

num_properties(model::LinearElastic) = num_properties(model.model) + 1
num_state_variables(model::LinearElastic) = 0

function property_names(model::LinearElastic)
    return ["density", property_names(model.model)...]
end

function cauchy_stress(
    model::LinearElastic,
    props, Z_old, Z_new, Δt, ε, θ
)
    props = module_props(model.model, props, 2)
    return cauchy_stress(model.model, props, ε, θ)
end

function helmholtz_free_energy(
    model::LinearElastic,
    props, Z_old, Z_new, Δt, ε, θ
)   props = module_props(model.model, props, 2)
    return helmholtz_free_energy(model.model, props, ε, θ)
end

function material_tangent(
    model::LinearElastic,
    props, Z_old, Z_new, Δt, ε, θ
)
    props = module_props(model.model, props, 2)
    return material_tangent(model.model, props, zero(Tensor{2, 3, eltype(ε), 9}), θ)
end

function spatial_tangent(
    model::LinearElastic,
    props, Z_old, Z_new, Δt, ε, θ
)
    props = module_props(model.model, props, 2)
    return spatial_tangent(model.model, props, ε, θ)
end

function p_wave_modulus(model::LinearElastic, props)
    props = module_props(model.model, props, 2)
    return p_wave_modulus(model.model, props)
end