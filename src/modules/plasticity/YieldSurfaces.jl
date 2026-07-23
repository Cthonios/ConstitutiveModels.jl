abstract type AbstractYieldSurface{I} <: AbstractConstitutiveModule end

# currently no models have kinematic hardening TODO
function has_kinematic_hardening(::AbstractYieldSurface)
    return false
end

# expected interface below
function effective_stress end

struct VonMises{I} <: AbstractYieldSurface{I}
    isotropic_hardening_model::I
end

function initialize_props(model::VonMises, inputs::Dict{String})
    return [
        initialize_props(model.isotropic_hardening_model, inputs)...
    ]
end

function initialize_state(model::VonMises)
    return [
        zero(SymmetricTensor{2, 3, Float64, 6}).data...,
        initialize_state(model.isotropic_hardening_model)...
    ]
end

num_properties(model::VonMises) = num_properties(model.isotropic_hardening_model)
num_state_variables(model::VonMises) = num_state_variables(model.isotropic_hardening_model)

function effective_stress(::VonMises, props, σ)
    return norm(dev(σ))
end

# # special case that doesn't require iteration
# function hardening_increment(
#     model::VonMises{I},
#     props_ih, μ, σ_eff, eqps_old::T, σ_b_old
# ) where {
#     I <: Union{LinearIsotropicHardening, NoIsotropicHardening},
#     T <: Number
# }
#     f = σ_eff - yield_surface_radius(model.isotropic_hardening_model, props_ih, eqps_old)
#     if f <= zero(T)
#         Δγ = zero(T)
#     else
#         K_prime = hardening_hessian(model.isotropic_hardening_model, props_ih, eqs_old)
#         Δγ = f / (T(2) * μ * (T(1) + K_prime / (T(3) * μ)))
#     end
#     return Δγ
# end

function hardening_increment(
    model::VonMises,
    props_ih, μ, σ_eff, eqps_old::T, σ_b_old
) where {
    T <: Number
}
    Δγ = zero(T)
    # eqps_new = eqps_old
    tol = 1e-8
    max_iter = 20
    
    # x = Δγ
    function f(x)
        eqps = eqps_old + sqrt(2. / 3.) * x
        K = hardening_gradient(model.isotropic_hardening_model, props_ih, eqps)
        g = -sqrt(2. / 3.) * K + σ_eff -
            (2μ * x) # TODO add kinematic hardening
        return g
    end

    function df(x)
        eqps = eqps_old + sqrt(2. / 3.) * x
        K_prime = hardening_hessian(model.isotropic_hardening_model, props_ih, eqps)
        dg = -2μ * (one(T) + (K_prime) / (3μ)) # TODO add kinematic hardening
        return dg
    end

    converged = false
    iter = 1
    while !converged
        g = f(Δγ)
        if abs(g)^2 < tol
            converged = true
            break
        end
        dg = df(Δγ)

        Δγ = Δγ - g / dg

        if iter >= max_iter
            @warn "MAx iteratations reached"
            break
        end
        iter += 1
    end
    return Δγ
end

function update_state(model::VonMises, props, μ::T, σ_tr, eqps_old, σ_b_old) where T <: Number
    σ_eff = effective_stress(model, props, σ_tr)
    props_ih = module_props(model.isotropic_hardening_model, props, 1)
    f_trial = σ_eff - yield_surface_radius(model.isotropic_hardening_model, props_ih, eqps_old)
    if f_trial <= zero(T)
        Δγ = zero(T)
        N = one(SymmetricTensor{2, 3, T, 6})
    else
        Δγ = hardening_increment(model, props_ih, μ, σ_eff, eqps_old, σ_b_old)
        N = dev(σ_tr) / norm(dev(σ_tr))
    end
    return Δγ, N
end
