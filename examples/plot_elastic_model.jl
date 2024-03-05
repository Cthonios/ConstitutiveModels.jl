using ConstitutiveModels
using Plots
using StaticArrays
using Tensors

inputs = Dict(
  "bulk modulus"  => 1000.0,#u"MPa",
  "shear modulus" => 1.0#u"MPa"
)
model = MechanicalModel(ConstitutiveModels.Hencky, inputs)

function run_loop!(Fs, σs, motion, model, inputs, λs)
  Δt = 0.0
  θ  = 0.0
  props = ConstitutiveModels.initialize_props(model, inputs)
  state_old = ConstitutiveModels.initialize_state(model)
  for λ in λs
    F = deformation_gradient(motion, λ, model, (props, Δt, θ, state_old))
    σ, state_new = cauchy_stress(model, props, Δt, F, θ, state_old)
    @show σ
    state_old = state_new
    push!(Fs, F)
    push!(σs, σ)
  end
end


motion = UniaxialStressDisplacementControl
λs = LinRange(1.0, 4., 100)

Fs = Tensor{2, 3, Float64, 9}[]
σs = Tensor{2, 3, Float64, 9}[]

@time run_loop!(Fs, σs, motion, model, inputs, λs)

F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1, 1], σs)

Ps = map((x, y) -> det(y) * dot(x, inv(y)'), σs, Fs)
P_11s = map(x -> x[1, 1], Ps)

p = plot(F_11s, σ_11s)
plot!(p, F_11s, P_11s)
display(P_11s)
# plot!(p, F_11s, σ_11_analytic)
