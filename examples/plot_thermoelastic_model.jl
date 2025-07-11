using ConstitutiveModels
using Plots

inputs = Dict(
  "Young's modulus"                  => 1.0,
  "Poisson's ratio"                  => 0.3,
  "coefficient of thermal expansion" => 0.001,
  "reference temperature"            => 60.0,
  "specific heat capacity"           => 1.0e-3,
  "thermal conductivity"             => 1.0e-4
)

motion = UniaxialStrain{Float64}()
λs = LinRange(1., 1.5, 101)

model = LinearThermoElastic()
props = initialize_props(model, inputs)
Δt = 0.0
θ = 90.0
Z = initialize_state(model)

@time ∇us, results = simulate(
  entropy, model, props, Δt, θ, Z, motion, λs
)
Fs = map(x -> x + one(x), ∇us)
F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1], results)
p = plot(F_11s, σ_11s)