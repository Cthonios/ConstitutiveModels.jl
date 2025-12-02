using ConstitutiveModels
using Enzyme
using Plots

inputs = Dict(
  "Young's modulus" => 1.0,
  "Poisson's ratio" => 0.49995,
  "n"               => 25
)

motion = UniaxialStressDisplacementControl{Float64}()
λs = LinRange(1.0, 8., 101)

model = NeoHookean()
# model = ArrudaBoyce()
props = initialize_props(model, inputs)
Δt = 0.0
θ  = 0.0
Z_old = initialize_state(model)
Z_new = initialize_state(model)

@time ∇us, σs, _ = simulate_material_point(
  pk1_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
)
Fs = map(x -> x + one(x), ∇us)
F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1, 1], σs)
p = plot(F_11s, σ_11s)
savefig(p, "elastic_plot.png")
