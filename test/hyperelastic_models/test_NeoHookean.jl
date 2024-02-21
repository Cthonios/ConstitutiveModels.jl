inputs_NeoHookean = Dict(
  "bulk modulus"  => 10.0,
  "shear modulus" => 1.0
)

@testset ExtendedTestSet "NeoHookean - Uniaxial stress, perfectly incompressible" begin
  model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
  @test props == [10.0, 1.0]
  λs = LinRange(0.5, 4., 100)
  ad_test(model, props, IsochoricUniaxialStress, λs)
end

@testset ExtendedTestSet "NeoHookean - Uniaxial strain" begin
  model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
  @test props == [10.0, 1.0]
  λs = LinRange(0.5, 4., 100)
  ad_test(model, props, UniaxialStrain, λs)
end

@testset ExtendedTestSet "NeoHookean - simple shear" begin
  model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
  @test props == [10.0, 1.0]
  λs = LinRange(-0.5, 0.5, 100)
  ad_test(model, props, SimpleShear, λs)
end

@testset ExtendedTestSet "NeoHookean - uniaxial stress, displacement controlled" begin
  model, props, state = MechanicalModel(NeoHookean, inputs_NeoHookean)
  @test props == [10.0, 1.0]
  λs = LinRange(0.5, 4.0, 100)
  ad_test(model, props, UniaxialStressDisplacementControl, λs)
end