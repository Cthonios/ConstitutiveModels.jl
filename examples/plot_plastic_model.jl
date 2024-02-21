using ConstitutiveModels
using Plots
using Tensors

props = Dict(
  "youngs modulus"            => 70e9,
  "poissons ratio"            => 0.3,
  "yield surface"             => "J2YieldSurface",
  "yield stress"              => 200.0e6,
  "isotropic hardening model" => "LinearIsotropicHardening",
  "hardening modulus"         => 200.0e6
  # "isotropic hardening model" => "SwiftVoceIsotropicHardening",
  # "hardening modulus"         => 200.0e6,
  # "A"                         => 200.0e6,
  # "n"                         => 20.0

)
model, props, state_old = MechanicalModel(LinearElastoPlasticity, props; type=SVector)

function run_loop!(Fs, σs, motion, model, props, state_init, states, λs)
  state_old = state_init
  state_new = copy(state_old)
  for λ in λs
    F = deformation_gradient(motion, model, props, state_old, λ, Tensor)
    # F = deformation_gradient(motion, λ)
    # σ = cauchy_stress(model, props, F)
    σ, props, state_new = ConstitutiveModels.cauchy_stress(model, props, F, state_old)

    state_old = state_new
    # @show σ
    push!(Fs, F)
    push!(σs, σ)
    push!(states, state_new)
  end
end


motion = UniaxialStressDisplacementControl
# motion = UniaxialStrain
λs = LinRange(1.0, 1.5, 1000)

Fs = Tensor{2, 3, Float64, 9}[]
σs = Tensor{2, 3, Float64, 9}[]
states = Vector{Float64}[]

@time run_loop!(Fs, σs, motion, model, props, state_old, states, λs)

F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1, 1], σs)
# display(σ_11s)
# Ps = map((x, y) -> det(y) * dot(x, inv(y)'), σs, Fs)
# P_11s = map(x -> x[1, 1], Ps)

p = plot(F_11s, σ_11s)
# plot!(p, F_11s, P_11s)
# plot!(p, F_11s, σ_11_analytic)
