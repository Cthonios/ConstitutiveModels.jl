struct ThermoHyperelastic{
    A <: AbstractHeatCapacity,
    B <: AbstractHeatConduction,
    C <: AbstractHyperelasticity
}
    heat_capacity::A
    heat_conduction::B
    hyperelasticity::C
end

function initialize_props(model::ThermoHyperelastic, inputs::PropertyInputs{T}) where T <: Number
    return [
        get_property(inputs, "density"),
        initialize_props(model.heat_capacity, inputs)...,
        initialize_props(model.heat_conduction, inputs)...,
        initialize_props(model.hyperelasticity, inputs)...
    ]
end
