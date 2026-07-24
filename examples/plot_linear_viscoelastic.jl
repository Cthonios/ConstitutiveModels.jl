using ConstitutiveModels
using Plots
using TabularFunctions

model = LinearViscoelastic(3)
props = Dict(
  "density"            => 1.0e3,
  "Young's modulus"    => 1.0e6,
  "Poisson's ratio"    => 0.3,
  "prony shear moduli" => [1.0e8, 2.0e8, 3.0e8],
  "relaxation times"   => [1.0, 10.0, 100.0]
)
λ_func = @piecewise_linear begin
    0.0, 1.00
    1.0, 1.05
    1000.0, 1.05
    1001.0, 1.10
    2001.0, 1.10
    2002.0, 1.00
end
motion = UniaxialStrain(λ_func)
out = simulate_material_point(cauchy_stress, model, props, motion, 2002.0)

t = map(x -> x.time, out)
ε = map(x -> linear_strain(motion, x)[1, 1], t)
σ = map(x -> x.material_output[1, 1], out)

λ = ε .+ 1

p1 = plot(λ, σ)
p2 = plot(t, σ)
# # plot!(p, λ, σ_an)
savefig(p1, "stress_strain.png")
savefig(p2, "stress_time.png")
