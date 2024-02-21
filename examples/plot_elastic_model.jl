using ConstitutiveModels
using Plots
using Tensors

props = Dict(
  "bulk modulus"  => 1000.0,#u"MPa",
  "shear modulus" => 1.0#u"MPa"
)
# model, props, state = NeoHookean(props)
model, props, state = MechanicalModel(NeoHookean, props)


# props = Dict(
#   "youngs modulus"            => 70e9,
#   "poissons ratio"            => 0.3,
#   "yield stress"              => 200.0e6,
#   "isotropic hardening model" => "LinearIsotropicHardening",
#   "hardening modulus"         => 200.0e6
#   # "isotropic hardening model" => "VoceIsotropicHardening",
#   # "A"                         => 200.0e6,
#   # "n"                         => 20

# )
# model, props, state = LinearElastoPlasticity(props)

function run_loop!(Fs, σs, motion, model, props, state, λs)
  for λ in λs
    F = deformation_gradient(motion, model, props, λ, Tensor)
    # P = pk1_stress(model, props, F)
    # σ, state = cauchy_stress(model, props, F, state)
    σ = cauchy_stress(model, props, F)
    # @show σ
    push!(Fs, F)
    push!(σs, σ)
  end
end


motion = UniaxialStressDisplacementControl
λs = LinRange(1.0, 4., 100)

Fs = Tensor{2, 3, Float64, 9}[]
σs = Tensor{2, 3, Float64, 9}[]

@time run_loop!(Fs, σs, motion, model, props, state, λs)

F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1, 1], σs)

Ps = map((x, y) -> det(y) * dot(x, inv(y)'), σs, Fs)
P_11s = map(x -> x[1, 1], Ps)

p = plot(F_11s, σ_11s)
plot!(p, F_11s, P_11s)
display(P_11s)
# plot!(p, F_11s, σ_11_analytic)
