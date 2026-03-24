struct Hyperelastic{
    NP, NS,
    StrainEnergyDensity <: AbstractHyperelastic{NP, NS}
} <: AbstractMechanicalModel{NP, NS}
    strain_energy_density::StrainEnergyDensity
end

function initialize_model(::Type{<:Hyperelastic}, inputs::Dict{String})
    strain_energy_density = eval(Symbol(inputs["strain energy density"]))
    return Hyperelastic(strain_energy_density())
end

function initialize_props(model::Hyperelastic, inputs::Dict{String})
    return initialize_props(model.strain_energy_density, inputs)
end

function cauchy_stress(
    model::Hyperelastic,
    props, Δt, Z_old, Z_new,
    ∇u, θ
)
    F = ∇u + one(∇u)
    return cauchy_stress(model.strain_energy_density, props, F, θ)
end

function helmholtz_free_energy(
    model::Hyperelastic,
    props, Δt, Z_old, Z_new,
    ∇u, θ
)
    F = ∇u + one(∇u)
    return strain_energy_density(model.strain_energy_density, props, F, θ)
end

function material_tangent(
    model::Hyperelastic,
    props, Δt, Z_old, Z_new,
    ∇u, θ
)
    F = ∇u + one(∇u)
    return material_tangent(model.strain_energy_density, props, F, θ)
end

function pk1_stress(
    model::Hyperelastic,
    props, Δt, Z_old, Z_new,
    ∇u, θ
)
    F = ∇u + one(∇u)
    return pk1_stress(model.strain_energy_density, props, F, θ)
end

p_wave_modulus(model::Hyperelastic, props) = p_wave_modulus(model.strain_energy_density, props)
