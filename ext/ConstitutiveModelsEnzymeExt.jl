module ConstitutiveModelsEnzymeExt

using ConstitutiveModels
using Enzyme
using StaticArrays
using Tensors

# Below implementation
# currently doesn't support the args...
# we'll need to modify the EnzymeAD init cache type

function _model_arg_types(::M, T=Float64) where {NP, NS, M <: ConstitutiveModels.AbstractConstitutiveModel{NP, NS}}
    return (
        Const{M}, 
        Const{SVector{NP, T}},
        Const{T},
        Active{Tensor{2, 3, T, 9}},
        Const{T},
        Const{SVector{NS, T}}
    )
end

function _model_tape_type(model::M, T=Float64) where {NP, NS, M <: ConstitutiveModels.AbstractConstitutiveModel{NP, NS}}
    arg_types = _model_arg_types(model, T)
    tape_type(
        ReverseSplitWithPrimal, 
        Const{typeof(ConstitutiveModels.helmholtz_free_energy)}, 
        Active,
        arg_types...
    )
end

function model_gradient(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props::SVector{NP, T},
    Δt::T,
    ∇u::Tensor{2, 3, T, 9},
    θ::T,
    Z::SVector{NS, T},
    args...
) where {NP, NS, T}
    arg_types = _model_arg_types(model, Float64)
    forward, reverse = autodiff_thunk(
        ReverseSplitWithPrimal, 
        Const{typeof(ConstitutiveModels.helmholtz_free_energy)},
        Active{Tuple{Float64, SVector{NS, Float64}}},
        arg_types...
    )
    tape, result, shadow_result = forward(
        Const(ConstitutiveModels.helmholtz_free_energy),
        Const(model), 
        Const(props),
        Const(Δt),
        Active(∇u),
        Const(θ),
        Const(Z)
    )
    out = reverse(
        Const(ConstitutiveModels.helmholtz_free_energy),
        Const(model), 
        Const(props),
        Const(Δt),
        Active(∇u),
        Const(θ),
        Const(Z),
        (1., ones(typeof(Z))),
        tape
    )
    return result, out
end

function model_gradient_deferred(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props::SVector{NP, T},
    Δt::T,
    ∇u::Tensor{2, 3, T, 9},
    θ::T,
    Z::SVector{NS, T},
    args...
) where {NP, NS, T}
    TapeType = _model_tape_type(model, Float64)
    arg_types = _model_arg_types(model, Float64)

    # first pass for gradients
    forward, reverse = autodiff_deferred_thunk(
        ReverseSplitWithPrimal, 
        TapeType, 
        Const{typeof(ConstitutiveModels.helmholtz_free_energy)},
        Active{Tuple{Float64, SVector{NS, Float64}}},
        arg_types...
    )

    tape, result, shadow_result = forward(
        Const(ConstitutiveModels.helmholtz_free_energy),
        Const(model), 
        Const(props),
        Const(Δt),
        Active(∇u),
        Const(θ),
        Const(Z)
    )
    out = reverse(
        Const(ConstitutiveModels.helmholtz_free_energy),
        Const(model), 
        Const(props),
        Const(Δt),
        Active(∇u),
        Const(θ),
        Const(Z),
        (1., ones(typeof(Z))),
        tape
    )
    # return out, result
    return out
end

function model_field_gradient_deferred(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props::SVector{NP, T},
    Δt::T,
    ∇u::Tensor{2, 3, T, 9},
    θ::T,
    Z::SVector{NS, T},
    args...
) where {NP, NS, T}
    return model_gradient_deferred(model, props, Δt, ∇u, θ, Z, args...)[1][4]
end

function model_field_hessian(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props::SVector{NP, T},
    Δt::T,
    ∇u::Tensor{2, 3, T, 9},
    θ::T,
    Z::SVector{NS, T},
    args...
) where {NP, NS, T}
    basis = ntuple(
        i -> Tensor{2, 3, Float64, 9}(
            ntuple(j -> ifelse(i == j, 1.0, 0.0), Val(9))
        ), 
        Val(9)
    )
    rows = autodiff(
        Forward,
        model_field_gradient_deferred,
        Const(model), 
        Const(props),
        Const(Δt),
        BatchDuplicated(∇u, basis),
        Const(θ),
        Const(Z)
    )[1]
    data = ntuple(i -> rows[(i-1) ÷ 9 + 1][(i-1) % 9 + 1], Val(81))
    return Tensor{4, 3, T, 81}(data)
end

function ConstitutiveModels.pk1_stress(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props::SVector{NP, T},
    Δt::T,
    ∇u::Tensor{2, 3, T, 9},
    θ::T,
    Z::SVector{NS, T},
    ::ConstitutiveModels.EnzymeAD,
    args...
) where {NP, NS, T}
    result, dresult = model_gradient(
        model,
        props,
        Δt,
        ∇u, θ, Z,
        args...
    )
    dresult[1][4], result[2]
end

function ConstitutiveModels.material_tangent(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props,
    Δt,
    ∇u, θ, Z,
    ::ConstitutiveModels.EnzymeAD,
    args...
) where {NP, NS}
    model_field_hessian(
        model,
        props,
        Δt,
        ∇u, θ, Z,
        # ad_cache,
        args...
    )
end

end # module
