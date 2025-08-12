using ConstitutiveModels
using Enzyme
using LaTeXStrings
using Plots
using Tensors

function solve()
    inputs = Dict(
        "bulk modulus"  => 100.,
        "shear modulus" => 1.,
        "n"             => 4
    )
    elastic_props = ConstitutiveModels.ElasticConstants(inputs)
    model = ConstitutiveModels.Hencky()
    # model = ArrudaBoyce()
    props = initialize_props(model, inputs)
    Z = initialize_state(model)
    θ = 0.
    Δt = 0.

    motion = UniaxialStrain()
    λs = LinRange(1., 1.5, 11)

    ∇us = Tensor{2, 3, Float64, 9}[]
    σs = SymmetricTensor{2, 3, Float64, 6}[]
    Zs = []
    σ_xx_ans = []
    σ_yy_ans = []
    for λ in λs
        ∇u = displacement_gradient(motion, λ)
        σ, Z = cauchy_stress(model, props, Δt, ∇u, θ, Z, ConstitutiveModels.EnzymeAD())
        # display(P)
        κ, μ = props
        
        σ_xx = (elastic_props.λ * log(λ) + 2 * μ * log(λ)) / λ
        σ_yy = (elastic_props.λ * log(λ)) / λ
        push!(∇us, ∇u)
        push!(σs, σ)
        push!(Zs, Z)
        push!(σ_xx_ans, σ_xx)
        push!(σ_yy_ans, σ_yy)
    end
    plot(motion, ∇us, σs, Zs)#, σ_xx_ans, σ_yy_ans)
end
solve()

# ∇u = one(Tensor{2, 3, Float64, 9})

# ConstitutiveModels.helmholtz_free_energy(
#     model, props, Δt, ∇u, θ, Z
# )

# @show P_enz = ConstitutiveModels.pk1_stress(
#     model, props, Δt, ∇u, θ, Z, 
#     ConstitutiveModels.EnzymeAD()
# )
# @show P_fd = ConstitutiveModels.pk1_stress(
#     model, props, Δt, ∇u, θ, Z, 
#     ConstitutiveModels.ForwardDiffAD()
# )[1]

# @show P_enz - P_fd
# # P

# @show A_enz = ConstitutiveModels.material_tangent(
#     model, props, Δt, ∇u, θ, Z, 
#     ConstitutiveModels.EnzymeAD()
# )
# @show A_fd = ConstitutiveModels.material_tangent(
#     model, props, Δt, ∇u, θ, Z, 
#     ConstitutiveModels.ForwardDiffAD()
# )[1]

# @show A_enz - A_fd

# x = rand(SVector{2, Float64})
# v = rand(SVector{2, Float64})

# function one_hot_basis(::Type{T}, ::Val{N}) where {T, N}
#     ntuple(i -> SVector{N, T}(ntuple(j -> ifelse(i == j, one(T), zero(T)), Val(N))), Val(N))
# end

# basis = one_hot_basis(Float64, Val(9))

# function f(x)
#     # return x[1] * x[1] + x[2] * x[1]
#     return norm(x)
# end

# function grad(x)
#     # autodiff_deferred(Reverse, Const(f), Active, Active(x))[1][1]
#     TapeType = tape_type(
#         ReverseSplitWithPrimal, 
#         Const{typeof(f)}, 
#         Active, 
#         # Duplicated{typeof(A)}, 
#         Active{typeof(x)}
#     )
#     forward, reverse = autodiff_deferred_thunk(
#         ReverseSplitWithPrimal, 
#         TapeType, 
#         Const{typeof(f)}, 
#         Active{Float64}, 
#         Active{typeof(x)}
#     )

#     tape, result, shadow_result  = forward(Const(f), Active(x))
#     ∂v = reverse(Const(f), Active(x), 1.0, tape)
#     ∂v
# end

# function hessian(x)
#     # basis = one_hot_basis(Float64, Val(2))
#     basis = ntuple(
#         i -> Tensor{2, 3, Float64, 9}(
#             ntuple(j -> ifelse(i == j, 1.0, 0.0), Val(9))
#         ), 
#         Val(9)
#     )
#     out = autodiff(
#         Forward,
#         grad,
#         BatchDuplicated(x, basis)
#     )
#     # dx
#     # SMatrix{4, 4, Float64, 16}
#     hcat(values(out[1])...)
# end

# x = rand(Tensor{2, 3, Float64, 9})
# grad_x = grad(x)
# hess_x = hessian(x)
