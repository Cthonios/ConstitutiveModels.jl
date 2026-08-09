using ConstitutiveModels
using Plots

model = FeFv()
props = Dict(
  "density"         => 1.0e3,
  "Young's modulus" => 1.0e6,
  "Poisson's ratio" => 0.4995,
  "G neq"           => 50.0e6,
  "relaxation time" => 100.0
)

p = plot()
for (rate, end_time) in zip([1.5e-1, 1.5e-2, 1.5e-3], [10.0, 100.0, 1000.0])
    λ_func = t -> 1 + rate * t
    motion = UniaxialStrain(λ_func)

    out = simulate_material_point(cauchy_stress, model, props, motion, end_time)

    t = map(x -> x.time, out)
    d = map(x -> x.kinematics[1, 1], out)
    σ = map(x -> x.material_output[1, 1], out)

    λ = d .+ 1

# p = plot(λ, σ)
    plot!(λ, σ)
end
savefig(p, "stress_strain.png")
