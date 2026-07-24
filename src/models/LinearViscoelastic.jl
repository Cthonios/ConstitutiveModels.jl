struct LinearViscoelastic{N} <: AbstractLinearElasticModel
    elastic_model::LinearElasticity
end

function LinearViscoelastic(num_prony_terms::Int)
    return LinearViscoelastic(Val{num_prony_terms}())
end

function LinearViscoelastic(::Val{NP}) where NP
    return LinearViscoelastic{NP}(
        LinearElasticity()
    )
end

function initialize_props(model::LinearViscoelastic, inputs::Dict{String})
    el = initialize_props(model.elastic_model, inputs)
    gs = inputs["prony shear moduli"]
    τs = inputs["relaxation times"]
    @assert length(gs) == num_prony_terms(model)
    @assert length(τs) == num_prony_terms(model)
    g0 = el[2] + reduce(+, gs)
    # el[2] = el[2] / g0
    γ∞ = el[2] / g0
    γs = gs / g0
    props = [inputs["density"], el..., γ∞]
    for n in 1:num_prony_terms(model)
        push!(props, γs[n])
        push!(props, τs[n])
    end
    return props
end

num_prony_terms(::LinearViscoelastic{N}) where N = N
num_properties(model::LinearViscoelastic) = 1 + num_properties(model.elastic_model) + 2 * num_prony_terms(model)
num_state_variables(model::LinearViscoelastic) = 6 * (num_prony_terms(model) + 1)

function pack_state!(state, ::LinearViscoelastic, n, h_n)
    range = 6 * (n - 1) + 1:6 * n
    state[range] .= h_n.data
end

# function state_variable_names(model::LinearViscoelastic)
#     return mapreduce(x -> state_variable_names(SymmetricTensor{2, 3, Float64, 6}, "h$x"), vcat, 1:num_prony_terms(model))
# end

function unpack_state(::LinearViscoelastic, state, n::Int)
    # range = start_index:start_index + 6 - 1
    range = 6 * (n - 1) + 1:6 * n
    h_n = SymmetricTensor{2, 3, eltype(state), 6}(@views state[range])
    return h_n
end

function cauchy_stress(
    model::LinearViscoelastic,
    props, Z_old, Z_new, Δt::T, ε, θ
) where T <: Number
    props_el = module_props(model.elastic_model, props, 2)
    σ0 = cauchy_stress(model.elastic_model, props_el, ε, θ)
    s0_old = unpack_state(model, Z_old, 1)
    s0_new = dev(σ0)
    pack_state!(Z_new, model, 1, s0_new)
    γ∞ = props[4]
    σ = vol(σ0) + γ∞ * s0_new
    for n in 1:num_prony_terms(model)
        γ, τ = props[4 + 2 * (n - 1) + 1], props[4 + 2 * (n - 1) + 2]
        h_old = unpack_state(model, Z_old, n + 1)
        h_new = exp(-Δt / τ) * h_old + exp(-Δt / 2τ) * (s0_new - s0_old)
        pack_state!(Z_new, model, n + 1, h_new)
        σ = σ + γ * h_new
    end
    return σ
end
