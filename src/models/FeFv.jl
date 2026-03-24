struct FeFv{
    NP, NS,
    EqStrainEnergyDensity <: AbstractHyperelastic
} <: AbstractMechanicalModel{NP, NS}
    eq_strain_energy_density::EqStrainEnergyDensity
    
end

# function FeFv(eq_strain_energy_density::AbstractHyperelastic)

# function initialize(model::FeFv)

function initialize_model(::Type{<:FeFv}, inputs::Dict{String})
    eq_model = eval(Symbol(inputs["equilibrium branch"]["strain energy density"]))()
    NP = num_properties(eq_model)
    NS = num_state_variables(eq_model)
    return FeFv{NP, NS, typeof(eq_model)}(eq_model)
end

function initialize_props(model::FeFv, inputs::Dict{String})
    eq_props = initialize_props(model.eq_strain_energy_density, inputs["equilibrium branch"])
    return vcat(eq_props)
end

function initialize_state(model::FeFv)
    return vcat(
        initialize_state(model.eq_strain_energy_density)
    )
end
