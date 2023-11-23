# ConstitutiveModels 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cthonios.github.io/ConstitutiveModels.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cthonios.github.io/ConstitutiveModels.jl/dev/) [![Build Status](https://github.com/Cthonios/ConstitutiveModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Cthonios/ConstitutiveModels.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/Cthonios/ConstitutiveModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Cthonios/ConstitutiveModels.jl)

ConstitutiveModels.jl aims to offer a general package for efficient implementation of general constitutive models with state.

Example LinearElastoPlasticity model

```
using ConstitutiveModels
using Plots
using Tensors

function run_loop!(Fs, σs, motion, model, props, state, λs)
  for λ in λs
    @show λ
    F = deformation_gradient(motion, model, props, state, λ)
    # P = pk1_stress(model, props, F)
    σ, state = cauchy_stress(model, props, F, state)
    push!(Fs, F)
    push!(σs, σ)
  end
end

props = Dict(
  "youngs modulus"            => 70e9,
  "poissons ratio"            => 0.3,
  "yield stress"              => 200.0e6,
  "isotropic hardening model" => "VoceIsotropicHardening",
  "A"                         => 200.0e6,
  "n"                         => 20

)
model, props, state = LinearElastoPlasticity(props)

motion = UniaxialStressDisplacementControl
λs = LinRange(1.0, 1.1, 100)

Fs = Tensor{2, 3, Float64, 9}[]
σs = Tensor{2, 3, Float64, 9}[]

run_loop!(Fs, σs, motion, model, props, state, λs)

F_11s = map(x -> x[1, 1], Fs)
σ_11s = map(x -> x[1, 1], σs)

p = plot(F_11s, σ_11s)

```