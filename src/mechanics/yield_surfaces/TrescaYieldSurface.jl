struct TrescaYieldSurface{NP, NS, IH} <: AbstractYieldSurface{NP, NS, IH}
    isotropic_hardening::IH
end

function TrescaYieldSurface()
    iso_hard = NoIsotropicHardening()
    return TrescaYieldSurface{1, 7, typeof(iso_hard)}(iso_hard)
end

function TrescaYieldSurface(
    iso_hard::AbstractIsotropicHardening
)
    NP = num_properties(iso_hard)
    NS = 7 + num_state_variables(iso_hard)
    return TrescaYieldSurface{NP, NS, typeof(iso_hard)}(iso_hard)
end

function initialize_props(model::TrescaYieldSurface, inputs::Dict{String})
    return initialize_props(model.isotropic_hardening, inputs)
end

function effective_stress(::TrescaYieldSurface, σ::SymmetricTensor{2, 3, T, 6}) where T <: Number
    σs = eigvals(σ)
    term1 = abs(σs[1] - σs[2])
    term2 = abs(σs[2] - σs[3])
    term3 = abs(σs[3] - σs[1])
    return 0.5 * max(term1, term2, term3)
end

function hardening_increment(
    model::TrescaYieldSurface,
    props, μ, σ_tr, α_old
)
    σ_eff = effective_stress(model, σ_tr)
    return hardening_increment(
        model.isotropic_hardening,
        props, μ, σ_eff, α_old
    )
end