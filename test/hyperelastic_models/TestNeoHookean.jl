inputs_NeoHookean = Dict(
  "density"               => 1.,
  "strain energy density" => "NeoHookean",
  "bulk modulus"          => 10.0,
  "shear modulus"         => 1.0
)

@testset ExtendedTestSet "NeoHookean - init" begin
  model = ConstitutiveModel(Hyperelastic, inputs_NeoHookean)
  props = initialize_properties(model, inputs_NeoHookean)
  state = initialize_state(model)

  @test length(props) == 3
  @test props[1] ≈ 1.
  @test props[2] ≈ 10.0
  @test props[3] ≈ 1.0 

  @test length(state) == 0
end

@testset ExtendedTestSet "NeoHookean - uniaxial strain" begin
  model = ConstitutiveModel(Hyperelastic, inputs_NeoHookean)
  props = initialize_properties(model, inputs_NeoHookean)

  K, μ = props[2], props[3]
  λs = LinRange(0.5, 4., 100)

  mat_states = MaterialState(Hyperelastic, inputs_NeoHookean, UniaxialStrain, λs)

  for (λ, mat_state) in zip(λs, mat_states)
    σ = mat_state.σ
    σ_xx = 0.5 * K * (λ - 1 / λ) + (2. / 3.) * (λ^2 - 1) * λ^(-5 / 3)
    σ_yy = 0.5 * K * (λ - 1 / λ) - (1. / 3.) * (λ^2 - 1) * λ^(-5 / 3)
    @test σ_xx ≈ σ[1, 1]
    @test σ_yy ≈ σ[2, 2]
    @test σ_yy ≈ σ[3, 3]
    @test 0.0 ≈ σ[1, 2]
    @test 0.0 ≈ σ[2, 3]
    @test 0.0 ≈ σ[3, 1]
  end
end

# @testset ExtendedTestSet "NeoHookean - Uniaxial stress, perfectly incompressible" begin
#   model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
#   @test props == [10.0, 1.0]
#   λs = LinRange(0.5, 4., 100)
#   ad_test(model, props, IsochoricUniaxialStress, λs)
# end

# @testset ExtendedTestSet "NeoHookean - Uniaxial strain" begin
#   model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
#   @test props == [10.0, 1.0]
#   λs = LinRange(0.5, 4., 100)
#   ad_test(model, props, UniaxialStrain, λs)
# end

# @testset ExtendedTestSet "NeoHookean - simple shear" begin
#   model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
#   @test props == [10.0, 1.0]
#   λs = LinRange(-0.5, 0.5, 100)
#   ad_test(model, props, SimpleShear, λs)
# end

# @testset ExtendedTestSet "NeoHookean - uniaxial stress, displacement controlled" begin
#   model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
#   @test props == [10.0, 1.0]
#   λs = LinRange(0.5, 4.0, 100)
#   ad_test(model, props, UniaxialStressDisplacementControl, λs)
# end