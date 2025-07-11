using ConstitutiveModels
using Plots
using Tensors

inputs = Dict(
  "Young's modulus"           => 70e9,
  "Poisson's ratio"           => 0.3,
  "yield surface"             => "VonMisesYieldSurface",
  "yield stress"              => 200.0e6,
  # "isotropic hardening model" => "LinearIsotropicHardening",
  "isotropic hardening model" => "LinearIsotropicHardening",
  "hardening modulus"         => 200.0e6
  # "isotropic hardening model" => "SwiftVoceIsotropicHardening",
  # "hardening modulus"         => 200.0e6,
  # "A"                         => 200.0e6,
  # "n"                         => 20.0

)

model = LinearElastoPlasticity(
  # VonMisesYieldSurface(
  TrescaYieldSurface(
    LinearIsotropicHardening()
  )
)
props = initialize_props(model, inputs)
Δt = 0.0
θ  = 0.0
Z = initialize_state(model)

motion = UniaxialStressDisplacementControl{Float64}()
λs = LinRange(1.0, 1.1, 101)

∇us, results = simulate(
  pk1_stress, model, props, Δt, θ, Z, motion, λs
)
Fs = map(x -> x + one(x), ∇us)
F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1][1, 1], results)

plot(F_11s, σ_11s)
