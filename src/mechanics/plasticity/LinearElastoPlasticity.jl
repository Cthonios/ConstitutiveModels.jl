struct LinearElastoPlasticity{NP, Yield} <: AbstractPlasticityModel{NP, 7, LinearElastic, Yield}
    elastic_model::LinearElastic
    yield_surface::Yield
end

function LinearElastoPlasticity(yield_surf_model::AbstractYieldSurface)
    elastic_model = LinearElastic()
    NP = num_properties(elastic_model) +
         num_properties(yield_surf_model)
    return LinearElastoPlasticity{NP, typeof(yield_surf_model)}(
        elastic_model, yield_surf_model
    )
end

function helmholtz_free_energy(
    model::LinearElastoPlasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    # unpack props
    elastic_props = elastic_properties(model, props)
    yield_props = yield_surface_properties(model, props)
    μ = elastic_props[2]

    # kinematics
    ε = linear_strain(∇u)

    # unpack state variables
    ε_p_old, α_old = unpack_state(model, Z_old)

    # calculate elastic trial stress
    ε_e_tr = ε - ε_p_old
    σ_e_tr = cauchy_stress(
        model.elastic_model, elastic_props, Δt,
        ε_e_tr, θ, Z_old, Z_new
    )

    # update state
    ε_p_new, α_new = update_state(model, yield_props, μ, σ_e_tr, ε_p_old, α_old)

    # update stuff to calculate energy
    ε_e = ε - ε_p_new

    # energies
    ψ_e = helmholtz_free_energy(
        model.elastic_model, 
        elastic_props, Δt,
        ε_e, θ, Z_old, Z_new
    )
    # TODO need to cleanup interface to hardening
    ψ_hard = energy(model.yield_surface.isotropic_hardening, yield_props, α_new)
    ψ = ψ_e + ψ_hard

    # pack state
    pack_state!(Z_new, model, ε_p_new, α_new)

    return ψ
end

function cauchy_stress(
    model::LinearElastoPlasticity,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    # unpack props
    elastic_props = elastic_properties(model, props)
    yield_props = yield_surface_properties(model, props)
    μ = elastic_props[2]

    # kinematics
    ε = linear_strain(∇u)

    # unpack state variables
    ε_p_old, α_old = unpack_state(model, Z_old)

    # calculate elastic trial stress
    ε_e_tr = ε - ε_p_old
    σ_e_tr = cauchy_stress(
        model.elastic_model, elastic_props, Δt,
        ε_e_tr, θ, Z_old, Z_new
    )

    # update state
    ε_p_new, α_new = update_state(model, yield_props, μ, σ_e_tr, ε_p_old, α_old)

    # update stuff to calculate stress
    ε_e = ε - ε_p_new

    # TODO need to fix this to be σ = σ_tr - 2μΔγN I think?
    # check Simo and Hughes
    σ_e = cauchy_stress(
        model.elastic_model, elastic_props, Δt,
        ε_e, θ, Z_old, Z_new
    )

    # pack state
    pack_state!(Z_new, model, ε_p_new, α_new)

    return σ_e
end

# function pk1_stress(
#     model::LinearElastoPlasticity,
#     props, Δt,
#     ∇u, θ, Z_old, Z_new,
#     ::ForwardDiffAD
# )
#     F = ∇u + one(typeof(∇u))
#     J = det(F)
#     σ = cauchy_stress(model, props, Δt, ∇u, θ, Z_old, Z_new)
#     return J * dot(σ, inv(F)')
# end

function pack_state!(state, ::LinearElastoPlasticity, ε_p, α)
    state[1:6] .= ForwardDiff.value.(ForwardDiff.value.(ε_p.data))
    state[7] = ForwardDiff.value(ForwardDiff.value(α))
    return nothing
end

function update_state(
    model::LinearElastoPlasticity, 
    yield_props, μ, σ_e_tr, ε_p_old, α_old
)
    # calculate hardening increment
    Δγ = hardening_increment(model.yield_surface, yield_props, μ, σ_e_tr, α_old)

    # radial return
    if Δγ > 0.0
        N = dev(σ_e_tr) / norm(dev(σ_e_tr))
        ε_p_new = ε_p_old + Δγ * N
        α_new   = α_old   + sqrt(2. / 3.) * Δγ
    else
        ε_p_new = ε_p_old
        α_new   = α_old
    end
    return ε_p_new, α_new
end

function unpack_state(::LinearElastoPlasticity, Z)
    # indices = SVector{6, Int}(1:6)
    ε_p = SymmetricTensor{2, 3, eltype(Z), 6}(@views Z[1:6])
    α = Z[7]
    return ε_p, α
end
