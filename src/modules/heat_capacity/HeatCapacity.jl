abstract type AbstractHeatCapacity <: AbstractConstitutiveModule end

# expected minimum interface
# function heat_capacity end
function helmholtz_free_energy end

struct SimpleHeatCapacity <: AbstractHeatCapacity
end

function initialize_props(::SimpleHeatCapacity, inputs::Dict{String})
    return [
        inputs["reference temperature"],
        inputs["specific heat capacity"]
    ]
end

num_properties(::SimpleHeatCapacity) = 2
num_state_variables(::SimpleHeatCapacity) = 0

function disspation(::SimpleHeatCapacity, props, ∇u, θ)
    return zero(typeof(θ))
end

function entropy(::SimpleHeatCapacity, props, ∇u, θ)
    θ_0, c = props[2], props[3]
    η = (c / θ_0) * (θ - θ_0)
    return η
end

function heat_capacity(::SimpleHeatCapacity, props, ∇u, θ)
    θ_0, c = props[1], props[2]
    return c * (θ / θ_0)
end

function helmholtz_free_energy(::SimpleHeatCapacity, props, ∇u, θ)
    θ_0, c = props[1], props[2]
    return -(c / (2. * θ_0)) * (θ - θ_0)^2
end
