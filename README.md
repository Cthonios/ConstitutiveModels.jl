# ConstitutiveModels 
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://cthonios.github.io/ConstitutiveModels.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://cthonios.github.io/ConstitutiveModels.jl/dev/) [![Build Status](https://github.com/Cthonios/ConstitutiveModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Cthonios/ConstitutiveModels.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/Cthonios/ConstitutiveModels.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Cthonios/ConstitutiveModels.jl)

ConstitutiveModels.jl aims to offer a general package for efficient implementation of general constitutive models with state.

Example ``Hyperelastic`` model

```
using ConstitutiveModels
using Plots

model = Hyperelastic(NeoHookean())
props = Dict(
  "density"         => 1.0e3,
  "Young's modulus" => 1.0e6,
  "Poisson's ratio" => 0.4995
)
λ_func = t -> 1 + 3t
motion = UniaxialStrain(λ_func)

out = simulate_material_point(cauchy_stress, model, props, motion, 1.0)

t = map(x -> x.time, out)
ε = map(x -> linear_strain(motion, x)[1, 1], t)
σ = map(x -> x.material_output[1, 1], out)

λ = ε .+ 1

p = plot(λ, σ)
savefig(p, "stress_strain.png")


```