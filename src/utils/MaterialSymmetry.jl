abstract type AbstractMaterialSymmetry{O} end

# expected interface
# function as_tensor(::AbstractMaterialSymmetry, props) end
# function initialize_props(::AbstractMaterialSymmetry, inputs, keys) end
# function num_properties(::AbstractMaterialSymmetry) = <int> end

# struct Cubic{O} <: AbstractMaterialSymmetry{O}
# end

# function initialize_props(::Cubic{4}, inputs, keys)

# end

struct Isotropy{O} <: AbstractMaterialSymmetry{O}
end

function as_tensor(::Isotropy{2}, props)
    return props[1] * one(SymmetricTensor{2, 3, eltype(props), 6})
end

function as_tensor(::Isotropy{4}, props)
    val_1, val_2 = props[1], props[2]
    δ(i, j) = i == j ? 1 : 0
    f = (i, j, k, l) -> val_1 * δ(i, j) * δ(k, l) + 
        val_2 * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
    return SymmetricTensor{4, 3, eltype(props), 36}(f)
end

function initialize_props(::Isotropy{2}, inputs, keys)
    return [inputs[keys[1]]]
end

function initialize_props(::Isotropy{4}, inputs, keys)
    return [inputs[keys[1]], inputs[keys[2]]]
end

num_properties(::Isotropy{2}) = 1
num_properties(::Isotropy{4}) = 2
