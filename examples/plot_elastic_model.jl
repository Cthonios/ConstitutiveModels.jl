using ConstitutiveModels
using Plots

model = Hyperelastic
props = Dict(
  "density"               => 2800.,
  "strain energy density" => "NeoHookean",
  "bulk modulus"          => 1000.0,
  "shear modulus"         => 1.0
)
motion = UniaxialStressDisplacementControl
λs = LinRange(1., 1.5, 500)

mat_states = MaterialState(model, props, motion, λs)

F_11s = map(x -> x.F[1, 1], mat_states)
P_11s = map(x -> x.P[1, 1] - x.P[2, 2], mat_states)

p = plot(F_11s, P_11s)
