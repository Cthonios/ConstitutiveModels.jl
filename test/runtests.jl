using Aqua
using ConstitutiveModels
using Tensors
using Test

function basis(k, l)
    return Tensor{2, 3, Float64, 9}((i, j) -> i == k && j == l ? 1 : 0)
end

"""
    fd_material_tangent(model, props, Δt, ∇u, θ, Z_old, Z_new; h=cbrt(eps(eltype(∇u))))

Compute a finite difference approximation of the first Piola-Kirchhoff stress
tensor **P** = ∂ψ/∂F by perturbing each component of the displacement gradient
∇u (i.e. the deformation gradient F = I + ∇u, or just F passed directly
depending on the model convention) with a central-difference stencil:

    A_ijkl ≈ [P(F + h·ei⊗ej⊗ek⊗el) - P(F - h·ei⊗ej⊗ek⊗el)] / (2h)
"""
function fd_material_tangent(
    model,
    props,
    Δt,
    Z_old,
    Z_new,
    ∇u::Tensor{2, 3, T, 9},
    θ,;
    h::T = cbrt(eps(T))
) where T
    data = ()
    for l in 1:3, k in 1:3
        δ = basis(k, l)

        ∇u_fwd = ∇u + h * δ
        ∇u_bwd = ∇u - h * δ

        P_fwd = pk1_stress(model, props, Δt, Z_old, Z_new, ∇u_fwd, θ)
        P_bwd = pk1_stress(model, props, Δt, Z_old, Z_new, ∇u_bwd, θ)

        ΔP = (P_fwd - P_bwd) / (2h)

        for j in 1:3, i in 1:3
            data = (data..., ΔP[i, j])
        end
    end
    return Tensor{4, 3, T, 81}(data)
end

"""
    fd_pk1_stress(model, props, Δt, ∇u, θ, Z_old, Z_new; h=cbrt(eps(eltype(∇u))))

Compute a finite difference approximation of the first Piola-Kirchhoff stress
tensor **P** = ∂ψ/∂F by perturbing each component of the displacement gradient
∇u (i.e. the deformation gradient F = I + ∇u, or just F passed directly
depending on the model convention) with a central-difference stencil:

    P_iJ ≈ [ψ(F + h·ei⊗ej) - ψ(F - h·ei⊗ej)] / (2h)
"""
function fd_pk1_stress(
    model,
    props,
    Δt,
    Z_old,
    Z_new,
    ∇u::Tensor{2, 3, T, 9},
    θ;
    h::T = cbrt(eps(T))
) where T
    data = ()
    for j in 1:3, i in 1:3
        δ = basis(i, j)

        # second order accurate
        ∇u_fwd = ∇u + h * δ
        ∇u_bwd = ∇u - h * δ
        ψ_fwd = helmholtz_free_energy(model, props, Δt, Z_old, Z_new, ∇u_fwd, θ)
        ψ_bwd = helmholtz_free_energy(model, props, Δt, Z_old, Z_new, ∇u_bwd, θ)
        val = (ψ_fwd - ψ_bwd) / (2h)
        # data = (data..., (ψ_fwd - ψ_bwd) / (2h))
        # fourth order accurate
        # ψ_fwd2 = helmholtz_free_energy(model, props, Δt, ∇u + 2h * δ, θ, Z_old, Z_new)
        # ψ_fwd1 = helmholtz_free_energy(model, props, Δt, ∇u +  h * δ, θ, Z_old, Z_new)
        # ψ_bwd1 = helmholtz_free_energy(model, props, Δt, ∇u -  h * δ, θ, Z_old, Z_new)
        # ψ_bwd2 = helmholtz_free_energy(model, props, Δt, ∇u - 2h * δ, θ, Z_old, Z_new)
        # val    = (-ψ_fwd2 + 8ψ_fwd1 - 8ψ_bwd1 + ψ_bwd2) / (12h)
        data = (data..., val)
    end

    return Tensors.Tensor{2, 3, T, 9}(data)
end

function test_strain_eq(::UniaxialStressDisplacementControl, εs, ε_xx_ans, ε_yy_ans)
    @assert length(εs) == length(ε_xx_ans)
    @assert length(εs) == length(ε_yy_ans)

    for (ε, ε_xx, ε_yy) in zip(εs, ε_xx_ans, ε_yy_ans)
        @test ε_xx ≈ ε[1, 1]
        @test ε_yy ≈ ε[2, 2]
        @test ε_yy ≈ ε[3, 3]
        @test 0.0 ≈ ε[1, 2]
        @test 0.0 ≈ ε[2, 3]
        @test 0.0 ≈ ε[3, 1]
    end
end

function test_stress_eq(::M, σs, σ_xx_ans, σ_yy_ans, σ_xy_ans; atol=1e-10, rtol=1e-10) where M <: Union{PureShearStrain, SimpleShear}
    @assert length(σs) == length(σ_xx_ans)
    @assert length(σs) == length(σ_yy_ans)
    @assert length(σs) == length(σ_xy_ans)
    for (σ, σ_xx, σ_yy, σ_xy) in zip(σs, σ_xx_ans, σ_yy_ans, σ_xy_ans)
        @test σ_xx ≈ σ[1, 1] atol=atol rtol=rtol
        @test σ_yy ≈ σ[2, 2] atol=atol rtol=rtol
        @test σ_yy ≈ σ[3, 3] atol=atol rtol=rtol
        @test σ_xy ≈ σ[1, 2] atol=atol rtol=rtol
        @test 0.0 ≈ σ[2, 3] atol=atol rtol=rtol
        @test 0.0 ≈ σ[3, 1] atol=atol rtol=rtol
    end
end

function test_stress_eq(::M, σs, σ_xx_ans, σ_yy_ans, σ_zz_ans, σ_xy_ans; atol=1e-10, rtol=1e-10) where M <: Union{PureShearStrain, SimpleShear}
    @assert length(σs) == length(σ_xx_ans)
    @assert length(σs) == length(σ_yy_ans)
    @assert length(σs) == length(σ_zz_ans)
    @assert length(σs) == length(σ_xy_ans)
    for (σ, σ_xx, σ_yy, σ_zz, σ_xy) in zip(σs, σ_xx_ans, σ_yy_ans, σ_zz_ans, σ_xy_ans)
        @test σ_xx ≈ σ[1, 1] atol=atol rtol=rtol
        @test σ_yy ≈ σ[2, 2] atol=atol rtol=rtol
        @test σ_zz ≈ σ[3, 3] atol=atol rtol=rtol
        @test σ_xy ≈ σ[1, 2] atol=atol rtol=rtol
        @test 0.0 ≈ σ[2, 3] atol=atol rtol=rtol
        @test 0.0 ≈ σ[3, 1] atol=atol rtol=rtol
    end
end

function test_stress_eq(::M, σs, σ_xx_ans, σ_yy_ans; atol=1e-10, rtol=1e-10) where M <: Union{UniaxialStrain, UniaxialStressDisplacementControl}
    @assert length(σs) == length(σ_xx_ans)
    @assert length(σs) == length(σ_yy_ans)
    for (σ, σ_xx, σ_yy) in zip(σs, σ_xx_ans, σ_yy_ans)
        @test σ_xx ≈ σ[1, 1] atol=atol rtol=rtol
        @test σ_yy ≈ σ[2, 2] atol=atol rtol=rtol
        @test σ_yy ≈ σ[3, 3] atol=atol rtol=rtol
        @test 0.0 ≈ σ[1, 2] atol=atol rtol=rtol
        @test 0.0 ≈ σ[2, 3] atol=atol rtol=rtol
        @test 0.0 ≈ σ[3, 1] atol=atol rtol=rtol
    end
end

@testset "Elastic constants" begin
    include("TestElasticConstants.jl")
end

@testset "Hyperelastic models" begin
    include("TestHyperelasticModels.jl")
end

@testset "Plasticity models" begin
    include("TestPlasticityModels.jl")
end

@testset "Utils" begin
    include("TestUtils.jl")
end

@testset "Aqua" begin
    Aqua.test_all(
        ConstitutiveModels,
        piracies=false
    )
end
