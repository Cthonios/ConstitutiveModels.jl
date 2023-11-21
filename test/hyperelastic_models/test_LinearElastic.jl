inputs_LinearElastic = Dict(
  "youngs modulus"  => 70.0e6,
  "poissons ratio" => 0.3
)

@testset ExtendedTestSet "LinearElastic - uniaxial strain unit test" begin
  model, props = LinearElastic(inputs_LinearElastic)
  λ_prop, μ = props[1], props[2]
  λs = LinRange(0.5, 1.5, 100)
  for λ in λs
    F = deformation_gradient(UniaxialStrain, λ)
    σ = cauchy_stress(model, props, F)
    ε_xx = λ - 1.
    σ_xx = λ_prop * ε_xx + 2. * μ * ε_xx
    σ_yy = λ_prop * ε_xx
    σ_zz = λ_prop * ε_xx
    @test σ_xx ≈ σ[1, 1]
    @test σ_yy ≈ σ[2, 2]
    @test σ_zz ≈ σ[3, 3]
    @test 0.0 ≈ σ[1, 2]
    @test 0.0 ≈ σ[2, 3]
    @test 0.0 ≈ σ[3, 1]
  end

  ad_test(model, props, UniaxialStrain, λs)
end

@testset ExtendedTestSet "LinearElastic - simple shear unit test" begin
  model, props = LinearElastic(inputs_LinearElastic)
  λ_prop, μ = props[1], props[2]
  λs = LinRange(-0.5, 0.5, 100)
  for λ in λs
    F = deformation_gradient(SimpleShear, λ)
    σ = cauchy_stress(model, props, F)
    ε_xy = λ
    @test 0.0 ≈ σ[1, 1]
    @test 0.0 ≈ σ[2, 2]
    @test 0.0 ≈ σ[3, 3]
    @test μ * ε_xy ≈ σ[1, 2]
    @test 0.0 ≈ σ[2, 3]
    @test 0.0 ≈ σ[3, 1]
  end

  ad_test(model, props, SimpleShear, λs)
end

# TODO add analytic solution test in below
@testset ExtendedTestSet "LinearElastic - uniaxial stress unit test" begin
  model, props = LinearElastic(inputs_LinearElastic)
  λs = LinRange(0.5, 1.5, 100)
  ad_test(model, props, UniaxialStressDisplacementControl, λs)
end