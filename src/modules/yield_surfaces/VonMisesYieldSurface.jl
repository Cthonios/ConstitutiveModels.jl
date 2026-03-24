struct VonMisesYieldSurface{NP, NS, IH} <: AbstractYieldSurface{NP, NS, IH}
    isotropic_hardening::IH
end

# default case for no hardening
function VonMisesYieldSurface()
    iso_hard = NoIsotropicHardening()
    return VonMisesYieldSurface{1, 7, typeof(iso_hard)}(iso_hard)
end

function VonMisesYieldSurface(
    iso_hard::AbstractIsotropicHardening
)
    NP = num_properties(iso_hard)
    NS = 7 + num_state_variables(iso_hard)
    return VonMisesYieldSurface{NP, NS, typeof(iso_hard)}(iso_hard)
end

function initialize_props(model::VonMisesYieldSurface, inputs::Dict{String})
    return initialize_props(model.isotropic_hardening, inputs)
end

function effective_stress(::VonMisesYieldSurface, σ::SymmetricTensor{2, 3, T, 6}) where T <: Number
    return norm(dev(σ))
end

function hardening_increment(
    model::VonMisesYieldSurface,
    props, μ, σ_tr, α_old
)
    σ_eff = effective_stress(model, σ_tr)
    return hardening_increment(
        model.isotropic_hardening,
        props, μ, σ_eff, α_old
    )
end
