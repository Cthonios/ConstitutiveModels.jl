module ConstitutiveModelsEnzymeExt

using ConstitutiveModels
using Enzyme
using StaticArrays
using Tensors

function _model_arg_types(::M, T=Float64) where {NP, NS, M <: ConstitutiveModels.AbstractConstitutiveModel{NP, NS}}
    return (
        Const{M}, 
        Const{SVector{NP, T}},
        Const{T},
        Active{Tensor{2, 3, T, 9}},
        Const{T},
        Const{Vector{T}},
        Const{Vector{T}}
    )
end

function model_gradient(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props::SVector{NP, T},
    Δt::T,
    ∇u::Tensor{2, 3, T, 9},
    θ::T,
    Z_old::AbstractArray{T, 1},
    Z_new::AbstractArray{T, 1},
    args...
) where {NP, NS, T}
    arg_types = _model_arg_types(model, Float64)
    forward, reverse = autodiff_thunk(
        ReverseSplitWithPrimal, 
        Const{typeof(ConstitutiveModels.helmholtz_free_energy)},
        # Active{Tuple{Float64, SVector{NS, Float64}}},
        Active{Tuple{Float64}},
        arg_types...
    )
    tape, result, shadow_result = forward(
        Const(ConstitutiveModels.helmholtz_free_energy),
        Const(model), 
        Const(props),
        Const(Δt),
        Active(∇u),
        Const(θ),
        Const(Z_old),
        Const(Z_new)
    )
    out = reverse(
        Const(ConstitutiveModels.helmholtz_free_energy),
        Const(model), 
        Const(props),
        Const(Δt),
        Active(∇u),
        Const(θ),
        Const(Z_old),
        Const(Z_new),
        # (1., ones(typeof(Z))),
        (1.,),
        tape
    )
    return result, out
    # autodiff(
    #     Reverse,
    #     Const(helmholtz_free_energy),
    #     Const(model),
    #     Const(props),
    #     Const(Δt),
    #     Active(∇u),
    #     Const(θ),
    #     Const(Z_old),
    #     Const(Z_new)
    # )
end

# function model_gradient_deferred(
#     model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
#     props::SVector{NP, T},
#     Δt::T,
#     ∇u::Tensor{2, 3, T, 9},
#     θ::T,
#     Z::SVector{NS, T},
#     args...
# ) where {NP, NS, T}
#     TapeType = _model_tape_type(model, Float64)
#     arg_types = _model_arg_types(model, Float64)

#     # first pass for gradients
#     forward, reverse = autodiff_deferred_thunk(
#         ReverseSplitWithPrimal, 
#         TapeType, 
#         Const{typeof(ConstitutiveModels.helmholtz_free_energy)},
#         Active{Tuple{Float64, SVector{NS, Float64}}},
#         arg_types...
#     )

#     tape, result, shadow_result = forward(
#         Const(ConstitutiveModels.helmholtz_free_energy),
#         Const(model), 
#         Const(props),
#         Const(Δt),
#         Active(∇u),
#         Const(θ),
#         Const(Z)
#     )
#     out = reverse(
#         Const(ConstitutiveModels.helmholtz_free_energy),
#         Const(model), 
#         Const(props),
#         Const(Δt),
#         Active(∇u),
#         Const(θ),
#         Const(Z),
#         (1., ones(typeof(Z))),
#         tape
#     )
#     # return out, result
#     return out
# end

# function model_field_gradient_deferred(
#     model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
#     props::SVector{NP, T},
#     Δt::T,
#     ∇u::Tensor{2, 3, T, 9},
#     θ::T,
#     Z::SVector{NS, T},
#     args...
# ) where {NP, NS, T}
#     return model_gradient_deferred(model, props, Δt, ∇u, θ, Z, args...)[1][4]
# end

# function model_field_hessian(
#     model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
#     props::SVector{NP, T},
#     Δt::T,
#     ∇u::Tensor{2, 3, T, 9},
#     θ::T,
#     Z::SVector{NS, T},
#     args...
# ) where {NP, NS, T}
#     basis = ntuple(
#         i -> Tensor{2, 3, Float64, 9}(
#             ntuple(j -> ifelse(i == j, 1.0, 0.0), Val(9))
#         ), 
#         Val(9)
#     )
#     rows = autodiff(
#         Forward,
#         model_field_gradient_deferred,
#         Const(model), 
#         Const(props),
#         Const(Δt),
#         BatchDuplicated(∇u, basis),
#         Const(θ),
#         Const(Z)
#     )[1]
#     data = ntuple(i -> rows[(i-1) ÷ 9 + 1][(i-1) % 9 + 1], Val(81))
#     return Tensor{4, 3, T, 81}(data)
# end

function ConstitutiveModels.pk1_stress(
    model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
    props::SVector{NP, T},
    Δt::T,
    ∇u::Tensor{2, 3, T, 9},
    θ::T,
    Z_old::AbstractArray{T, 1},
    Z_new::AbstractArray{T, 1},
    ::ConstitutiveModels.EnzymeAD,
    args...
) where {NP, NS, T}
    result, dresult = model_gradient(
        model,
        props,
        Δt,
        ∇u, θ, Z_old, Z_new,
        args...
    )
    # dresult[1][4], result[2]
end

# function ConstitutiveModels.material_tangent(
#     model::ConstitutiveModels.AbstractMechanicalModel{NP, NS},
#     props,
#     Δt,
#     ∇u, θ, Z,
#     ::ConstitutiveModels.EnzymeAD,
#     args...
# ) where {NP, NS}
#     model_field_hessian(
#         model,
#         props,
#         Δt,
#         ∇u, θ, Z,
#         # ad_cache,
#         args...
#     )
# end

end # module
