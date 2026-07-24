using ConstitutiveModels
using Plots

model = Hyperelastic(NeoHookean())
props = Dict(
  "density"         => 1.0e3,
  "Young's modulus" => 1.0e6,
  "Poisson's ratio" => 0.4995
)
λ_func = t -> 1 + 3t
motion = UniaxialStrain(λ_func)

out = simulate_material_point(cauchy_stress, model, props, motion, 1.0)

t = map(x -> x.time, out)
ε = map(x -> linear_strain(motion, x)[1, 1], t)
σ = map(x -> x.material_output[1, 1], out)

λ = ε .+ 1

p = plot(λ, σ)
savefig(p, "stress_strain.png")
