struct LinearElastoplastic{I} <: AbstractLinearElasticModel
    elastic_model::LinearElasticity
    yield_surface::VonMises{I}
end

function LinearElastoplastic(isotropic_hardening)
    return LinearElastoplastic(
        LinearElasticity(),
        VonMises(isotropic_hardening)
    )
end

function initialize_props(model::LinearElastoplastic, inputs::Dict{String})
    return [
        inputs["density"],
        initialize_props(model.elastic_model, inputs)...,
        initialize_props(model.yield_surface, inputs)...
    ]
end

num_properties(model::LinearElastoplastic) = num_properties(model.elastic_model) + num_properties(model.yield_surface) + 1
num_state_variables(model::LinearElastoplastic) = 6 + num_state_variables(model.yield_surface)

function pack_state!(Z, model::LinearElastoplastic, ε_p, eqps, σ_b)
    # @show ε_p
    Z[1:6] .= ForwardDiff.value.(ε_p.data)
    Z[7] = ForwardDiff.value(eqps)
    if has_kinematic_hardening(model.yield_surface)
        @assert false "Implement me"
    end
    return nothing
end

function state_variable_names(::LinearElastoplastic)
    return [
        state_variable_names(SymmetricTensor{2, 3, T, 6}, "ep")...,
        "eqps"
    ]
end

function unpack_state(model::LinearElastoplastic, Z)
    p_indices = SVector{6, Int}(1:6)
    ε_p = SymmetricTensor{2, 3, eltype(Z), 6}(@views Z[p_indices])
    eqps = Z[7]
    if has_kinematic_hardening(model.yield_surface)
        @assert false "Implement me"
    else
        σ_b = zero(SymmetricTensor{2, 3, eltype(Z), 6})
    end
    return ε_p, eqps, σ_b
end

# function helmholtz_free_energy(
#     model::LinearElastoplastic,
#     props, Z_old, Z_new, Δt, ε, θ
# )
#     props_el = module_props(model.elastic_model, props, 2)
#     props_pl = module_props(model.yield_surface, props, 4)
#     ε_p_old, eqps_old, σ_b_old = unpack_state(model, Z_old)
#     ε_e_tr = ε - ε_p_old
#     σ_tr = cauchy_stress(model.elastic_model, props_el, ε_e_tr, θ)
#     ε_p_new, eqps_new, σ_b_new = update_state(model.yield_surface, props_pl, μ, σ_tr, eqps_old, σ_b_old)
#     ε_e_new = ε - ε_p_new
#     ψ_e = helmholtz_free_energy(model.elastic_model, props_el, ε_e_new, θ)
#     # ψ_h = hardening_energy(model.yield_surface.isotropic_hardening_model, )

#     pack_state!(Z_new, model, ε_p_new, eqps_new, σ_b_new)

#     ψ = ψ_e + ψ_h
#     return ψ
# end

function cauchy_stress(
    model::LinearElastoplastic,
    props, Z_old, Z_new, Δt::T, ε, θ
) where T <: Number
    props_el = module_props(model.elastic_model, props, 2)
    μ = props_el[2]
    props_pl = module_props(model.yield_surface, props, 4)
    ε_p_old, eqps_old, σ_b_old = unpack_state(model, Z_old)
    ε_e_tr = ε - ε_p_old
    σ_tr = cauchy_stress(model.elastic_model, props_el, ε_e_tr, θ)
    Δγ, N = update_state(model.yield_surface, props_pl, μ, σ_tr, eqps_old, σ_b_old)
    ε_p_new = ε_p_old + Δγ * N
    eqps_new = eqps_old + sqrt(T(2) / T(3)) * Δγ
    # TODO implement back stress
    σ_b_new = zero(SymmetricTensor{2, 3, T, 6})
    pack_state!(Z_new, model, ε_p_new, eqps_new, σ_b_new)
    σ = σ_tr - 2μ * Δγ * N
    return σ
end

function spatial_tangent(
    model::LinearElastoplastic,
    props, Z_old, Z_new, Δt::T, ε, θ
) where T <: Number
    props_el = module_props(model.elastic_model, props, 2)
    μ = props_el[2]
    props_pl = module_props(model.yield_surface, props, 4)
    ε_p_old, eqps_old, σ_b_old = unpack_state(model, Z_old)
    ε_e_tr = ε - ε_p_old
    σ_tr = cauchy_stress(model.elastic_model, props_el, ε_e_tr, θ)
    Δγ, N = update_state(model.yield_surface, props_pl, μ, σ_tr, eqps_old, σ_b_old)
    ε_p_new = ε_p_old + Δγ * N
    eqps_new = eqps_old + sqrt(T(2) / T(3)) * Δγ
    # TODO implement back stress
    σ_b_new = zero(SymmetricTensor{2, 3, T, 6})
    pack_state!(Z_new, model, ε_p_new, eqps_new, σ_b_new)
    # σ = σ_tr - 2μ * Δγ * N
    # return σ
    C = spatial_tangent(model.elastic_model, props_el, ε, θ)

    # TODO make below call yield surface methods that wrap hardening methods
    K_prime = hardening_hessian(model.yield_surface.isotropic_hardening_model, props_pl, eqps_new)
    H_prime = zero(T) # TODO
    term_1 = -2μ * Δγ / norm(dev(σ_tr))
    term_2 = one(T) / (one(T) + (K_prime + H_prime) / (3μ)) + term_1
    I = one(SymmetricTensor{2, 3, T, 6})
    II = one(SymmetricTensor{4, 3, T, 36})
    C = C + 2μ * term_1 * (II - (one(T) / T(3)) * otimes(I, I)) - 2μ * term_2 * otimes(N, N)
    return C
end

function p_wave_modulus(::LinearElastoplastic, props) 
    props = module_props(model.elastic_model, props, 2)
    return p_wave_modulus(model.elastic_model, props)
end