abstract type AbstractConstitutiveModel{NP, NS} end
num_properties(::AbstractConstitutiveModel{NP, NS}) where {NP, NS} = NP
num_state_variables(::AbstractConstitutiveModel{NP, NS}) where {NP, NS} = NS

"""
All models are assumed to take in a type
tag for the model and a Dict of props

This method we likely can't statically compile
due to the Any
"""
function initialize_props(
    model, 
    inputs::Dict{String}
)
    new_inputs = Dict{Symbol, Any}()
    for (key, val) in inputs
      new_inputs[Symbol(key)] = val
    end
    return initialize_props(model, new_inputs)
end

"""
Default state constructor to just return zeros
"""
function initialize_state(
    ::AbstractConstitutiveModel{NP, NS},
    float_type = Float64
) where {NP, NS}
    return zero(SVector{NS, float_type})
end

@inline function unpack_props(
    ::AbstractConstitutiveModel{NP1, NS},
    props::SVector{NP2, T},
    start_index::Int
) where {NP1, NP2, NS, T <: Number}
    indices = SVector{NP1, typeof(NP1)}(start_index:start_index + NP1 - 1)
    return props[indices]
end

# function unpack_state(
#     ::AbstractConstitutiveModel{NP, NS1},
#     state::SVector{NS2, T},
#     start_index::Int
# ) where {NP, NS1, NS2, T <: Number}
#     indices = SVector{NS1, typeof(NS1)}(start_index:start_index + NS1 - 1)
#     return state[indices]
# end

# other then these methods...
# nothing really makes sense in general
# e.g. diffusion doesn't really have
# a "helmholtz_free_energy"
# so other interfaces will be defined
# for different physics types

# however, the following convention will be followed
# for a physics type with a method called for example
# ```physics_qoi```, we will have the following method
# signature convention
# ```physics_qoi(model, props, Δt, ∇u_1, ∇u_2, ...,  θ, Z, args...)```
# where model is the instance of a given ```AbstractConstitutiveModel```
# ```props``` is a pre-initialized flat array of properties
# ```Δt``` is the time step, ```∇u_1```, etc. 
# are the field gradients (this convention will break for hessian theories),
# ```θ``` is the absolute temperature, and ```Z``` is the
# collection of current, e.g. old, state variables.
# 
# each of these methods is assumed to return a Tensor
# e.g. rank 0, 1, 2, 4, etc., and a new set of updated
# state variables
#
# Even for models without state variables, e.g. elasticity,
# this is still necessary

# Additional things for AD
abstract type AbstractADType end
struct ForwardDiffAD <: AbstractADType 
end
struct EnzymeAD <: AbstractADType
end

# abstract type AbstractProperty end

const Properties{N, T} = SVector{N, T} where {N, T}
