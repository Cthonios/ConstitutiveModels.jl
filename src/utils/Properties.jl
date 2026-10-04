function print_property_keys(io::IO, inputs)
    for key in keys(inputs)
        println(io, "  ", key)
    end
end

# function print_property_keys(io::IO, inputs::PropertyInputs{T}) where T <: Number
#     for key in keys(inputs)
#         println(io, "  ", key)
#     end
# end

const PropertyType{T} = Union{T, String} where T <: Number

@inline _fill!(d::Dict{String, Any}) = d
@inline function _fill!(
    d::Dict{String, Any}, p::Pair{String}, rest::Vararg{Pair{String}, N}
) where N
    d[p.first] = p.second
    return _fill!(d, rest...)
end

struct PropertyInputs{T <: Number}
    data::Dict{String, Any}

    function PropertyInputs()
        return PropertyInputs{Float64}()
    end

    function PropertyInputs{T}() where T <: Number
        data = Dict{String, Any}()
        new{T}(data)
    end

    function PropertyInputs(data::Dict{String})
        return PropertyInputs{Float64}(data)
    end

    function PropertyInputs{T}(data::Dict{String}) where T <: Number
        new{T}(data)
    end

    function PropertyInputs(
        first::Pair{String}, pairs::Vararg{Pair{String}, N}
    ) where N
        return PropertyInputs{Float64}(first, pairs...)
    end

    function PropertyInputs{T}(
        first::Pair{String}, pairs::Vararg{Pair{String}, N}
    ) where {T <: Number, N}
        return new{T}(_fill!(Dict{String, Any}(), first, pairs...))
    end
end

function Base.get(inputs::PropertyInputs{T}, key::String, default) where T <: Number
    return get(inputs.data, key, default)
end

function Base.getindex(inputs::PropertyInputs{T}, key::String) where T <: Number
    return inputs.data[key]
end

function Base.haskey(inputs::PropertyInputs{T}, key::String) where T <: Number
    return haskey(inputs.data, key)
end

function Base.keys(inputs::PropertyInputs{T}) where T <: Number
    return keys(inputs.data)
end

function Base.setindex!(inputs::PropertyInputs{T}, val, key::String) where T <: Number
    inputs.data[key] = val
    return nothing
end

function Base.values(inputs::PropertyInputs{T}) where T <: Number
    return values(inputs.data)
end

function Base.show(io::IO, inputs::PropertyInputs{T}) where T <: Number
    for key in keys(inputs)
        println(io, "  ", key)
    end
    # println(io, inputs.data)
end

struct BadNestedPropertyError{T <: Number} <: Exception
    inputs::PropertyInputs{T}
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

struct PropertyNotFoundError{T <: Number} <: Exception
    inputs::PropertyInputs{T}
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
        return val::Float64
    elseif isa(val, Dict)
        # one of many potential cases
        if haskey(val, "mean") && haskey(val, "variance")
            dist = Normal(val["mean"], val["variance"])
            return rand(dist)::Float64
        else
            bad_nested_property_input(val, key)
        end
    else
        # default case
        return default
    end
end

function get_property(inputs::PropertyInputs{T}, key::String, type::Type{String}, default = nothing) where T <: Number
    if default === nothing
        if !haskey(inputs, key)
            property_not_found_error(inputs, key)
        end
        return inputs[key]::String
    else
        return get(inputs, key, default)
    end
end

function get_property(inputs::PropertyInputs{T}, key::String, type::Type{T}, default = nothing) where T <: Number
    if default === nothing
        if !haskey(inputs, key)
            property_not_found_error(inputs, key)
        end
        return inputs[key]::T
    else
        return get(inputs, key, default)
    end
end

function get_property(inputs::PropertyInputs{T}, key::String, type::Type{<:Dict{String}}) where T <: Number
    # val = get(inputs, key, default)
    # if val === nothing
    #     property_not_found_error(inputs, key)
    # else
    #     return val
    # end
    # return val
    return inputs.data[key]
end
