using ConstitutiveModels
using Plots
using Tensors

inputs = Dict(
  "youngs modulus"            => 70e9,
  "poissons ratio"            => 0.3,
  "yield surface"             => "J2YieldSurface",
  "yield stress"              => 200.0e6,
  # "isotropic hardening model" => "LinearIsotropicHardening",
  "isotropic hardening model" => "LinearIsotropicHardening",
  "hardening modulus"         => 200.0e6
  # "isotropic hardening model" => "SwiftVoceIsotropicHardening",
  # "hardening modulus"         => 200.0e6,
  # "A"                         => 200.0e6,
  # "n"                         => 20.0

)
model = MechanicalModel(LinearElastoPlasticity, inputs)

function run_loop!(Fs, σs, motion, model, inputs, λs)
  Δt = 0.0
  θ  = 0.0
  props = ConstitutiveModels.initialize_props(model, inputs)
  state_old = ConstitutiveModels.initialize_state(model)

  for λ in λs
    @show λ
    F = deformation_gradient(motion, λ, model, (props, Δt, θ, state_old))
    σ, state_new = ConstitutiveModels.cauchy_stress(model, props, Δt, F, θ, state_old)
    state_old = state_new
    push!(Fs, F)
    push!(σs, σ)
    # push!(states, state_new)
  end
end


motion = UniaxialStressDisplacementControl
λs = LinRange(1.0, 1.05, 1000)

Fs = Tensor{2, 3, Float64, 9}[]
σs = Tensor{2, 3, Float64, 9}[]
states = Vector{Float64}[]

@time run_loop!(Fs, σs, motion, model, inputs, λs)

F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1, 1], σs)
# display(σ_11s)
# Ps = map((x, y) -> det(y) * dot(x, inv(y)'), σs, Fs)
# P_11s = map(x -> x[1, 1], Ps)

p = plot(F_11s, σ_11s)
# plot!(p, F_11s, P_11s)
# plot!(p, F_11s, σ_11_analytic)
