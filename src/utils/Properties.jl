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

function get_property(inputs::Dict{String}, key::String, default = nothing)
    val = get(inputs, key, default)
    if val === nothing
        property_not_found_error(inputs, key)
    elseif isa(val, Number)
        # scalar value case
        # if isa(val, Float64)
        #     return val::Float64
        # return val::Float64
        return val
    elseif isa(val, Dict)
        # one of many potential cases
        if haskey(val, "mean") && haskey(val, "variance")
            dist = Normal(val["mean"], val["variance"])
            return rand(dist)
        else
            bad_nested_property_input(val, key)
        end
    else
        # default case
        return default
    end
end
