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

# function LinearElastoPlasticity(inputs::Dict{String})
#     yield_surf_type = eval(Meta.parse(inputs["yield surface"]))

#     elastic_model = LinearElastic()
#     yield_surf_model = yield_surf_type()
#     n_props = num_properties(elastic_model) +
#               num_properties(yield_surf_model)
#     return LinearElastoPlasticity{n_props, typeof(yield_surf_model)}(
#         elastic_model,
#         yield_surf_model
#     )
# end

function helmholtz_free_energy(
    model::LinearElastoPlasticity,
    props, Δt,
    ∇u, θ, Z
)

    # unpack props
    elastic_props = elastic_properties(model, props)
    yield_props = yield_surface_properties(model, props)
    μ = elastic_props[2]

    # kinematics
    ε = symmetric(∇u)

    # unpack state variables
    ε_p_old, α_old = unpack_state(model, Z)

    # calculate elastic trial stress
    ε_e_tr = ε - ε_p_old
    σ_e_tr, _ = cauchy_stress(
        model.elastic_model, elastic_props, Δt,
        ε_e_tr, θ, SVector{0, typeof(θ)}()
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

    # update stuff to calculate energy
    ε_e = ε - ε_p_new

    # energies
    ψ_e, _ = helmholtz_free_energy(
        model.elastic_model, 
        elastic_props, Δt,
        ε_e, θ, SVector{0, typeof(θ)}()
    )
    # TODO need to cleanup interface to hardening
    ψ_hard = energy(model.yield_surface.isotropic_hardening, yield_props, α_new)
    ψ = ψ_e + ψ_hard

    # pack state
    Z = pack_state(model, ε_p_new, α_new)

    return ψ, Z
end

function pack_state(::LinearElastoPlasticity, ε_p, α)
    return SVector{7, eltype(ε_p)}(ε_p.data..., α)
end

function unpack_state(::LinearElastoPlasticity, Z)
    indices = SVector{6, Int}(1:6)
    ε_p = SymmetricTensor{2, 3, eltype(Z), 6}(Z[indices])
    α = Z[7]
    return ε_p, α
end
