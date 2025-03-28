inputs_LinearElastic = Dict(
  "density"               => 1.0, 
  "strain energy density" => "LinearElastic",
  "youngs modulus"        => 70.0e9,
  "poissons ratio"        => 0.3
)

# @testset ExtendedTestSet "LinearElastic - uniaxial strain unit test" begin
  model = ConstitutiveModel(Hyperelastic, inputs_LinearElastic)
  props = initialize_properties(model, inputs_LinearElastic)
  # state_old = initialize_state(model)

  # mat_states = Mater 
  # Δt = 0.0
  # θ = 0.0
  λ_prop, μ = props[2], props[3]
  λs = LinRange(0.5, 1.5, 100)

  mat_states = MaterialState(Hyperelastic, inputs_LinearElastic, UniaxialStrain, λs)

  for (λ, mat_state) in zip(λs, mat_states)
    # F = ConstitutiveModels.deformation_gradient(UniaxialStrain, λ)
    # σ, state_new = cauchy_stress(model, props, Δt, F, θ, state_old)
    σ = mat_state.σ
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

  # ad_test(model, props, state, UniaxialStrain, λs)
# end

# @testset ExtendedTestSet "LinearElastic - simple shear unit test" begin
#   model = MechanicalModel(LinearElastic, inputs_LinearElastic)
#   props = ConstitutiveModels.initialize_props(model, inputs_LinearElastic)
#   state_old = ConstitutiveModels.initialize_state(model)
#   Δt = 0.0
#   θ = 0.0
#   λ_prop, μ = props[1], props[2]
#   λs = LinRange(-0.5, 0.5, 100)
#   for λ in λs
#     F = ConstitutiveModels.deformation_gradient(SimpleShear, λ)
#     σ, state_new = cauchy_stress(model, props, Δt, F, θ, state_old)
#     ε_xy = λ
#     @test 0.0 ≈ σ[1, 1]
#     @test 0.0 ≈ σ[2, 2]
#     @test 0.0 ≈ σ[3, 3]
#     @test μ * ε_xy ≈ σ[1, 2]
#     @test 0.0 ≈ σ[2, 3]
#     @test 0.0 ≈ σ[3, 1]
#   end

#   # ad_test(model, props, state, SimpleShear, λs)
# end

# TODO add analytic solution test in below
# @testset ExtendedTestSet "LinearElastic - uniaxial stress unit test" begin
#   model, props, state = LinearElastic(inputs_LinearElastic)
#   λs = LinRange(0.5, 1.5, 100)
#   ad_test(model, props, state, UniaxialStressDisplacementControl, λs)
# end