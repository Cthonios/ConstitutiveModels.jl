struct Hypoelastic{
    H <: AbstractObjectiveStressRate
} <: AbstractHypoelasticModel
    rate::H
end

function initialize_props(::Hypoelastic, inputs::Dict{String})
    el = ElasticConstants(inputs)
    return [get_property(inputs, "density"), el.λ, el.μ]
end

function initialize_state(model::Hypoelastic)
    state = zeros(num_state_variables(model))
    state[7] = 1.0
    return state
end

num_properties(::Hypoelastic) = 3
num_state_variables(::Hypoelastic) = 7

function pack_state!(Z, ::Hypoelastic, σ, J)
    Z[1:6] .= σ.data
    Z[7] = J
    return nothing
end

function unpack_state(::Hypoelastic, state)
    σ = SymmetricTensor{2, 3, Float64, 6}(@views state[1:6])
    J = state[7]
    return σ, J
end

function cauchy_stress(
    model::Hypoelastic,
    props, Z_old, Z_new, Δt, L, θ
)
    # props = module_props(model.rate, 2)
    props = SVector{2, eltype(props)}(props[2], props[3])
    # return cauchy_stress(model.rate, props, Δt, σ_old, L, θ)
    σ, J = unpack_state(model, state)
    σ, J = update_state(model.rate, props, Δt, J, σ, L, θ)
end
