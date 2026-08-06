struct ElasticConstantsError <: Exception
    msg::String
end

function Base.showerror(io::IO, e::ElasticConstantsError)
    names = [
        "bulk modulus",
        "Lamé's first constant",
        "Poisson's ratio",
        "shear modulus",
        "Young's Modulus"
    ]
    println(io, typeof(e), "\n", e.msg)
    println(io, "Possible values include:")
    for name in names
        println(io, "  $name")
    end
end

function elastic_constants_error(msg::String)
    throw(ElasticConstantsError(msg))
end

"""
$(TYPEDEF)
$(TYPEDFIELDS)
"""
struct ElasticConstants{T <: Number}
    E::T
    ν::T
    κ::T
    λ::T
    μ::T
end

# TODO add unit init
"""
Allowable parameters names are the following

- bulk modulus
- Lamé's first constant
- Poisson's ratio
- shear modulus
- Young's modulus
$(TYPEDSIGNATURES)
"""
function ElasticConstants(params::Dict{String})
    if haskey(params, "Young's modulus") == true
        E = get_property(params, "Young's modulus")
        if haskey(params, "Poisson's ratio") == true
            ν = get_property(params, "Poisson's ratio")
            κ = E / 3(1 - 2ν)
            λ = E * ν / (1 + ν) / (1 - 2ν)
            μ = E / 2(1 + ν)
        elseif haskey(params, "bulk modulus") == true
            κ = get_property(params, "bulk modulus")
            ν = (3κ - E) / 6κ
            λ = (3κ * (3κ - E)) / (9κ - E)
            μ = 3κ * E / (9κ - E)
        elseif haskey(params, "Lamé's first constant") == true
            λ = get_property(params, "Lamé's first constant")
            R = sqrt(E^2 + 9λ^2 + 2E * λ)
            ν = 2λ / (E + λ + R)
            κ = (E + 3λ + R) / 6
            μ = (E - 3λ + R) / 4
        elseif haskey(params, "shear modulus") == true
            μ = get_property(params, "shear modulus")
            ν = E / 2μ - 1
            κ = E * μ / 3(3μ - E)
            λ = μ * (E - 2μ) / (3μ - E)
        else
            elastic_constants_error("Two elastic constants are required but only elastic modulus found")
        end
    elseif haskey(params, "Poisson's ratio") == true
        ν = get_property(params, "Poisson's ratio")
        if haskey(params, "bulk modulus") == true
            κ = get_property(params, "bulk modulus")
            E = 3κ * (1 - 2ν)
            λ = 3κ * ν / (1 + ν)
            μ = 3κ * (1 - 2ν) / 2(1 + ν)
        elseif haskey(params, "Lamé's first constant") == true
            λ = get_property(params, "Lamé's first constant")
            E = λ * (1 + ν) * (1 - 2ν) / ν
            κ = λ * (1 + ν) / 3ν
            μ = λ * (1 - 2ν) / 2ν
        elseif haskey(params, "shear modulus") == true
            μ = get_property(params, "shear modulus")
            E = 2μ * (1 + ν)
            κ = 2μ * (1 + ν) / 3(1 - 2ν)
            λ = 2μ * ν / (1 - 2ν)
        else
            elastic_constants_error("Two elastic constants are required but only Poisson's ratio found")
        end
    elseif haskey(params, "bulk modulus") == true
        κ = params["bulk modulus"]::Float64
        if haskey(params, "Lamé's first constant") == true
            λ = params["Lamé's first constant"]::Float64
            E = 9κ * (κ - λ) / (3κ - λ)
            ν = λ / (3κ - λ)
            μ = 3(κ - λ) / 2
        elseif haskey(params, "shear modulus") == true
            μ = params["shear modulus"]::Float64
            E = 9κ * μ / (3κ + μ)
            ν = (3κ - 2μ) / 2(3κ + μ)
            λ = κ - 2μ / 3
        else
            elastic_constants_error("Two elastic constants are required but only bulk modulus found")
        end
    elseif haskey(params, "Lamé's first constant") == true
        λ = params["Lamé's first constant"]::Float64
        if haskey(params, "shear modulus") == true
            μ = params["shear modulus"]::Float64
            E = μ * (3λ + 2μ) / (λ + μ)
            ν = λ / 2(λ + μ)
            κ = λ + 2μ / 3
        else
            elastic_constants_error("Two elastic constants are required but only Lamé's first constant found")
        end
    elseif haskey(params, "shear modulus") == true
        elastic_constants_error("Two elastic constants are required but only shear modulus found")
    else
        elastic_constants_error("Two elastic constants are required but none found")
    end
    return ElasticConstants(E, ν, κ, λ, μ)
end

function Base.show(io::IO, e::ElasticConstants)
    println(io, "Elastic constants:")
    println(io, "  bulk modulus          = $(e.κ)")
    println(io, "  Lamé's first constant = $(e.λ)")
    println(io, "  Poisson's ratio       = $(e.ν)")
    println(io, "  shear modulus         = $(e.μ)")
    println(io, "  Young's Modulus       = $(e.E)")
end

function to_dict(e::ElasticConstants{T}) where T <: Number
    d = Dict{String, T}(
        "bulk modulus"          => e.κ,
        "Lamé's first constant" => e.λ,
        "Poisson's ratio"       => e.ν,
        "shear modulus"         => e.μ,
        "Young's modulus"       => e.E
    )
    return d
end
