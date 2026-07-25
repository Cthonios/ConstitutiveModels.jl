abstract type AbstractThermalExpansion <: AbstractConstitutiveModule end
num_state_variables(::AbstractThermalExpansion) = 0

struct LinearThermalExpansion{S <: AbstractMaterialSymmetry} <: AbstractThermalExpansion
    symmetry::S
end

function initialize_props(model::LinearThermalExpansion, inputs::Dict{String})
    θ_0 = inputs["reference temperature"]
    # kind of hacky below for now but works
    if isa(model.symmetry, Isotropy{2})
        elastic_props = ElasticConstants(inputs)
        α = get_property(inputs, "coefficient of thermal expansion")
        β = -3 * elastic_props.κ * α
        return [θ_0, β]
    else
        @assert false "Currently unsupported material symmetry type $(model.symmetry)"
    end
end

num_properties(model::LinearThermalExpansion) = 1 + num_properties(model.symmetry)

function cauchy_stress(model::LinearThermalExpansion, props, ::SymmetricTensor{2, 3, T, 6}, θ) where T <: Number
    θ_0 = props[1]
    M_props = module_props(model.symmetry, props, 2)
    M = as_tensor(model.symmetry, M_props)
    return (θ - θ_0) * M
end

function entropy(model::LinearThermalExpansion, props, ε::SymmetricTensor{2, 3, T, 6}, θ) where T <: Number
    θ_0 = props[1]
    M_props = module_props(model.symmetry, props, 2)
    M = as_tensor(model.symmetry, M_props)
    return -θ_0 * dcontract(M, ε)
end

function helmholtz_free_energy(model::LinearThermalExpansion, props, ε::SymmetricTensor{2, 3, T, 6}, θ) where T <: Number
    θ_0 = props[1]
    M_props = module_props(model.symmetry, props, 2)
    M = as_tensor(model.symmetry, M_props)
    return (θ - θ_0) * dcontract(M, ε)
end
