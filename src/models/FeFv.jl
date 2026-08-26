# struct FeFv <: AbstractHyperelasticModel
#     model_eq::Hencky
# end

# function FeFv()
#     return FeFv(Hencky())
# end

# function initialize_props(model::FeFv, inputs::Dict{String})
#     return [
#         get_property(inputs, "density"),
#         initialize_props(model.model_eq, inputs)...,
#         get_property(inputs, "G neq"),
#         get_property(inputs, "relaxation time")
#     ]
# end

# function initialize_state(::FeFv)
#     return collect(one(Tensor{2, 3, Float64, 9}).data)
# end

# num_properties(model::FeFv) = 1 + num_properties(model.model_eq) + 2
# num_state_variables(model::FeFv) = 9

# # function pack_state!(Z, ::FeFv, Fv)
# #     Z[1:9] .= ForwardDiff.value.(Fv.data)
# #     return nothing
# # end
# function pack_state!(Z, ::FeFv, Fv)
#     Z[1:9] .= ForwardDiff.value.(ForwardDiff.value.(Fv.data))
#     return nothing
# end

# function property_names(model::FeFv)
#     return ["density", property_names(model.model), "G neq", "relaxation time"]
# end

# function state_variable_names(::FeFv)
#     return state_variable_names(Tensor{2, 3, Float64, 9}, "Fv")
# end

# function unpack_state(::FeFv, Z)
#     p_indices = SVector{9, Int}(1:9)
#     Fv = Tensor{2, 3, eltype(Z), 9}(@views Z[p_indices])
#     return Fv
# end

# function helmholtz_free_energy(
#     model::FeFv,
#     props, Z_old, Z_new, Δt, ∇u, θ
# )
#     # equilibrium response
#     props_eq = module_props(model.model_eq, props, 2)
#     ψ_eq = helmholtz_free_energy(model.model_eq, props_eq, ∇u, θ)

#     # calculate trial kinematics
#     F = one(∇u) + ∇u
#     Fv_old = unpack_state(model, Z_old)
#     Fe_trial = dot(F, inv(Fv_old))
#     Ce_trial = tdot(Fe_trial)
#     Ee_trial = 0.5 * log(Ce_trial)

#     # get state increment
#     G_neq, τ = props[4], props[5]
#     integration_factor = one(Δt) / (one(Δt) + Δt / τ)
#     # Ee_dev = dcontract(Ee_trial, Ee_trial)
#     Ee_dev = dev(Ee_trial)
#     ΔE_v = Δt * integration_factor * Ee_dev / τ
#     Fv_new = dot(exp(Δt * ΔE_v), Fv_old)
#     pack_state!(Z_new, model, Fv_new)

#     # neq strain energy
#     Ee_new = Ee_trial - ΔE_v
#     Ee_dev = dev(Ee_new)
#     ψ_neq = G_neq * dcontract(Ee_dev, Ee_dev)

#     # dissipation potential
#     Dv = ΔE_v / Δt
#     η = G_neq * τ
#     ϕ = η * dcontract(Dv, Dv)

#     return ψ_eq + ψ_neq + Δt * ϕ
# end

# # function pk1_stress(
# #     model::FeFv,
# #     props, Z_old, Z_new, Δt, ∇u, θ
# # )
# #     # equilibrium response
# #     props_eq = module_props(model.model_eq, props, 2)
# #     P_eq = pk1_stress(model.model_eq, props_eq, ∇u, θ)
    
# #     # calculate trial kinematics
# #     F = one(∇u) + ∇u
# #     Fv_old = unpack_state(model, Z_old)
# #     Fe_trial = dot(F, inv(Fv_old))
# #     Ce_trial = tdot(Fe_trial)
# #     Ee_trial = 0.5 * log(Ce_trial)
    
# #     # get state increment
# #     G_neq, τ = props[4], props[5]
# #     integration_factor = one(Δt) / (one(Δt) + Δt / τ)
# #     # Ee_dev = dcontract(Ee_trial, Ee_trial)
# #     Ee_dev = dev(Ee_trial)
# #     ΔE_v = Δt * integration_factor * Ee_dev / τ
# #     Fv_new = dot(exp(Δt * ΔE_v), Fv_old)
# #     pack_state!(Z_new, model, Fv_new)

# #     # neq response
# #     Ee_new = Ee_trial - ΔE_v
# #     Ee_dev = dev(Ee_new)
# #     # ψ_neq = G_neq * Ee_new
    


# #     return P_eq
# # end

# # function material_tangent(
# #     model::FeFv,
# #     props, Z_old, Z_new, Δt, ∇u, θ
# # )
# #     props_eq = module_props(model.model_eq, props, 2)
# #     A_eq = material_tangent(model.model_eq, props_eq, ∇u, θ)
# #     # TODO do visco stuff
# #     return A_eq
# # end

# function p_wave_modulus(model::FeFv, props)
#     props_eq = module_props(model.model_eq, props, 2)
#     return p_wave_modulus(model.model_eq, props_eq)
# end

struct FeFv <: AbstractHyperelasticModel
    model_eq::Hencky
end

function FeFv()
    return FeFv(Hencky())
end

function initialize_props(model::FeFv, inputs::Dict{String})
    return [
        get_property(inputs, "density"),
        initialize_props(model.model_eq, inputs)...,
        get_property(inputs, "G neq"),
        get_property(inputs, "relaxation time")
    ]
end

function initialize_state(::FeFv)
    return collect(one(Tensor{2, 3, Float64, 9}).data)
end

num_properties(model::FeFv) = 1 + num_properties(model.model_eq) + 2
num_state_variables(model::FeFv) = 9

function pack_state!(Z, ::FeFv, Fv)
    Z[1:9] .= ForwardDiff.value.(Fv.data)
    return nothing
end

function property_names(model::FeFv)
    return ["density", property_names(model.model_eq)..., "G neq", "relaxation time"]
end

function state_variable_names(::FeFv)
    return state_variable_names(Tensor{2, 3, Float64, 9}, "Fv")
end

function unpack_state(::FeFv, Z)
    p_indices = SVector{9, Int}(1:9)
    Fv = Tensor{2, 3, eltype(Z), 9}(@views Z[p_indices])
    return Fv
end

function neq_trial_state(model::FeFv, props, Z_old, Δt, ∇u)
    F = one(∇u) + ∇u
    Fv_old = unpack_state(model, Z_old)
    Fe_trial = dot(F, inv(Fv_old))
    Ce_trial = tdot(Fe_trial)
    Ee_trial = 0.5 * log(Ce_trial)
    Ee_dev = dev(Ee_trial)

    G_neq, τ = props[4], props[5]
    β = Δt / (τ + Δt)                     # = Δt * integration_factor / τ
    ΔE_v = β * Ee_dev
    η = G_neq * τ
    C = G_neq * (1 - β)^2 + η * β^2 / Δt  # effective (zero-bulk) shear modulus

    return (Fv_old = Fv_old, Fe_trial = Fe_trial, Ee_dev = Ee_dev, ΔE_v = ΔE_v, C = C)
end

function effective_shear_props(::Hencky, C)
    return [zero(C), C]
end

function helmholtz_free_energy(
    model::FeFv,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    # equilibrium response
    props_eq = module_props(model.model_eq, props, 2)
    ψ_eq = helmholtz_free_energy(model.model_eq, props_eq, ∇u, θ)

    # non-equilibrium (Maxwell) response
    st = neq_trial_state(model, props, Z_old, Δt, ∇u)
    Fv_new = dot(exp(st.ΔE_v), st.Fv_old)
    pack_state!(Z_new, model, Fv_new)
    ψ_neq_total = st.C * dcontract(st.Ee_dev, st.Ee_dev)

    return ψ_eq + ψ_neq_total
end

function pk1_stress(
    model::FeFv,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    props_eq = module_props(model.model_eq, props, 2)
    P_eq = pk1_stress(model.model_eq, props_eq, ∇u, θ)

    st = neq_trial_state(model, props, Z_old, Δt, ∇u)
    Fv_new = dot(exp(st.ΔE_v), st.Fv_old)
    pack_state!(Z_new, model, Fv_new)

    B = transpose(inv(st.Fv_old))
    ∇u_e = st.Fe_trial - one(st.Fe_trial)
    props_neq_eff = effective_shear_props(model.model_eq, st.C)
    Σ = pk1_stress(model.model_eq, props_neq_eff, ∇u_e, θ)
    P_neq = dot(Σ, B)

    return P_eq + P_neq
end

function material_tangent(
    model::FeFv,
    props, Z_old, Z_new, Δt, ∇u, θ
)
    props_eq = module_props(model.model_eq, props, 2)
    A_eq = material_tangent(model.model_eq, props_eq, ∇u, θ)

    st = neq_trial_state(model, props, Z_old, Δt, ∇u)
    B = transpose(inv(st.Fv_old))
    ∇u_e = st.Fe_trial - one(st.Fe_trial)
    props_neq_eff = effective_shear_props(model.model_eq, st.C)
    𝔸 = material_tangent(model.model_eq, props_neq_eff, ∇u_e, θ)

    A_neq = Tensor{4, 3}((i, L, j, M) -> sum(
        𝔸[i, N, j, Q] * B[N, L] * B[Q, M]
        for N in 1:3, Q in 1:3
    ))

    return A_eq + A_neq
end

function p_wave_modulus(model::FeFv, props)
    props_eq = module_props(model.model_eq, props, 2)
    return p_wave_modulus(model.model_eq, props_eq)
end
