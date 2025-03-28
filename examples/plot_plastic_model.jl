using ConstitutiveModels
using Plots


E = 70e9
ν = 0.3
G = E / (2 * (1 + ν))
K = E / (3 * (1 - 2 * ν))

model = LinearElastoPlastic
# model = FeFp
props = Dict(
  "density"             => 2800.,
  "youngs modulus"      => 70e9,
  "poissons ratio"      => 0.3,
  "bulk modulus"        => K,
  "shear modulus"       => G,
  "yield surface"       => "VonMises",
  "yield stress"        => 200.0e6,
  "isotropic hardening" => "LinearIsotropicHardening",
  "hardening modulus"   => 200.0e6
)
motion = UniaxialStressDisplacementControl
# motion = UniaxialStrain
λs = LinRange(1., 1.05, 500)

mat_states = MaterialState(model, props, motion, λs)

F_11s = map(x -> x.F[1, 1], mat_states)
P_11s = map(x -> x.P[1, 1] - x.P[2, 2], mat_states)

p = plot(F_11s, P_11s)
