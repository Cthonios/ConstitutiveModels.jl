using Aqua
using ConstitutiveModels
using Tensors
using Test

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

function test_stress_eq(::SimpleShear, σs, σ_xx_ans, σ_yy_ans, σ_xy_ans; atol=1e-10, rtol=1e-10)
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

@testset "Aqua" begin
    Aqua.test_all(
        ConstitutiveModels,
        piracies=false
    )
end
