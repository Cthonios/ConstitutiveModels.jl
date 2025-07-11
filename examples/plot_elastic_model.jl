using ConstitutiveModels
using Plots

inputs = Dict(
  "Young's modulus" => 1.0,
  "Poisson's ratio" => 0.49995
)

motion = UniaxialStressDisplacementControl{Float64}()
λs = LinRange(1.0, 4., 101)

model = NeoHookean()
props = initialize_props(model, inputs)
Δt = 0.0
θ  = 0.0
Z = initialize_state(model)

@time ∇us, results = simulate(
  pk1_stress, model, props, Δt, θ, Z, motion, λs
)
Fs = map(x -> x + one(x), ∇us)
F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1][1, 1], results)
p = plot(F_11s, σ_11s)
