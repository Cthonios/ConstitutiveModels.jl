function print_property_keys(io::IO, inputs::Dict{String})
    for key in keys(inputs)
        println(io, "  ", key)
    end
end

struct BadNestedPropertyError{D <: Dict} <: Exception
    inputs::D
    key::String
end

function Base.showerror(io::IO, err::BadNestedPropertyError)
    println(io, "Bad nested property input for key $(err.key).")
    println(io, "Nested input that threw the error is $(err.inputs).")
    println(io, "Currently supported nested input formats include:\n")
    println(io, "1. Normal distribution samples")
    println(io, "  \"key\" => Dict(\n    \"mean\"     => <real>\n    \"variance\" => <real>\n  )")
end

function bad_nested_property_input(inputs, key)
    throw(BadNestedPropertyError(inputs, key))
end

struct PropertyNotFoundError{D <: Dict} <: Exception
    inputs::D
    key::String
end

function Base.showerror(io::IO, err::PropertyNotFoundError)
    println(io, "Could not find property in inputs with key $(err.key).")
    println(io, "Properties passed by user are:")
    print_property_keys(io, err.inputs)
end

function property_not_found_error(inputs, key)
    throw(PropertyNotFoundError(inputs, key))
end

function _nested_helper(inputs::Dict{String}, key::String, default, ::Type{T}) where T
    if !haskey(inputs, key)
        return default
    end

    val = inputs[key]

    if val isa T
        return val::T
    elseif val isa Dict
        nested = val::Dict{String, T}
        if !haskey(nested, "mean") ||
           !haskey(nested, "variance")
            bad_nested_property_input(inputs, key)
        end

        mean = nested["mean"]
        variance = nested["variance"]

        if !(mean isa T) || !(variance isa T)
            bad_nested_property_input(inputs, key)
        end

        μ = mean::T
        σ² = variance::T

        return rand(Normal(μ, sqrt(σ²)))
    else
        bad_nested_property_input(inputs, key)
    end
end

function get_complex_property(inputs::Dict{String}, key::String, default::ComplexF64)
    return _nested_helper(inputs, key, default, ComplexF64)
end

function get_complex_property(inputs::Dict{String}, key::String)
    if !haskey(inputs, key)
        property_not_found_error(inputs, key)
    end
    return get_complex_property(inputs, key, 0.0 + 0.0im)
end

function get_real_property(inputs::Dict{String}, key::String, default::Float64)
    return _nested_helper(inputs, key, default, Float64)
end

function get_real_property(inputs::Dict{String}, key::String)
    if !haskey(inputs, key)
        property_not_found_error(inputs, key)
    end
    return get_real_property(inputs, key, 0.0)
end

function get_string_property(inputs::Dict{String}, key::String, default = nothing)
    val = get(inputs, key, default)
    if val isa String
        return val
    else
        property_not_found_error(inputs, key)
    end
end
