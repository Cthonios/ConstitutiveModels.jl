using ConstitutiveModels
using Plots

model = LinearElastoplastic(Voce())
props = Dict(
  "density"                          => 1.0e3,
  "Young's modulus"                  => 70e9,
  "Poisson's ratio"                  => 0.3,
  "hardening modulus"                => 0e6,
  "hardening exponent"               => 20,
  "saturation stress"                => 200e6,
  "yield stress"                     => 200e6,
  #
  "coefficient of thermal expansion" => 1e-3,
  "reference temperature"            => 30.0,
  "specific heat capacity"           => 1000.0,
  "thermal conductivity"             => 1e-4
)
λ_func = t -> 1 + 0.05t
# motion = UniaxialStrain(λ_func)
motion = UniaxialStressDisplacementControl(λ_func)
out = simulate_material_point(cauchy_stress, model, props, motion, 1.0)

ε = map(x -> x.kinematics[1, 1], out)
σ = map(x -> x.material_output[1, 1], out)

λ = ε .+ 1

p = plot(λ, σ)
savefig(p, "stress_strain.png")


