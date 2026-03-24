# useful for setting up different types of symmetry
# for different tensors
abstract type AbstractMaterialSymmetry end
struct Isotropic <: AbstractMaterialSymmetry
end

function initialize_props(
    ::Type{<:SymmetricTensor{2, 3}}, ::Type{<:Isotropic},
    inputs::Dict{String}, key
)
    val = inputs[key]
    return [val]
end

function initialize_props(
    ::Type{<:SymmetricTensor{4, 3}}, ::Type{<:Isotropic},
    inputs::Dict{String}, key_1, key_2
)
    val_1 = inputs[key]
    val_2 = inputs[key]
    return [val_1, val_2]
end

function initialize_tensor(
    type::Type{<:SymmetricTensor{2, 3}}, ::Type{<:Isotropic},
    props
)
    val = props[1]
    return val * one(type)
end

function initialize_tensor(
    type::Type{<:SymmetricTensor{4, 3}}, ::Type{<:Isotropic},
    props
)
    val_1 = props[1]
    val_2 = props[2]
    f = (i,j,k,l) -> val_1*δ(i,j)*δ(k,l) + val_2*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    return type(f)
end

function initialize_tensor(
    type::Type{<:SymmetricTensor{4, 3}}, ::Type{<:Isotropic},
    elastic_constants::ElasticConstants
)
    val_1 = elastic_constants.λ
    val_2 = elastic_constants.μ
    f = (i,j,k,l) -> val_1*δ(i,j)*δ(k,l) + val_2*(δ(i,k)*δ(j,l) + δ(i,l)*δ(j,k))
    return type(f)
end
