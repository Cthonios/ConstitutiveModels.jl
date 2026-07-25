# NeoHookean
```@autodocs
Modules = [ConstitutiveModels]
Order   = [:type, :function]
Pages   = ["NeoHookean.jl"]
```

## Pure Shear Strain
### Analytic Solution
``\mathbf{\sigma}_{11} = \frac{\mu}{3}\left[\frac{1}{2}\left(\lambda^2 + \lambda^{-2}\right) - 1\right]``

``\mathbf{\sigma}_{22} = \mathbf{\sigma}_{11}``

``\mathbf{\sigma}_{33} = \frac{\mu}{3}\left(2 - \lambda^2 + \lambda^{-2}\right)``

``\mathbf{\sigma}_{12} = \frac{\mu}{2}\left(\lambda^2 - \lambda^{-2}\right)``

All other components are zero

### Verification
Here is a comparison of an analytic solution to the uniaxial stress boundary value problem in displacement control.
```@example
using ConstitutiveModels
using Plots

function neohookean_pure_shear_strain()
    inputs = Dict(
        "density"         => 1.0,
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = Hyperelastic(NeoHookean())
    motion = PureShearStrain(t -> 1 + 1.25t)
    out = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)

    props = initialize_props(model, inputs)
    ts = map(x -> x.time, out)
    ∇us = map(x -> x.kinematics, out)
    λs = map(t -> 1 + 1.25t, ts)
    σs = map(x -> x.material_output, out)
    Zs = map(x -> x.state, out)

    μ = props[3]
    σ_11s_an = (μ / 3) * (0.5 * (λs.^2 .+ λs.^(-2)) .- 1)
    σ_33s_an = (μ / 3) * (2 .- λs.^2 .- λs.^(-2))
    σ_12s_an = (μ / 2) * (λs.^2 .- λs.^(-2))

    plot(motion, ∇us, σs, Zs, σ_11s_an, σ_33s_an, σ_12s_an)
end
neohookean_pure_shear_strain()
```

## Simple Shear
### Analytic Solution
``\mathbf{\sigma}_{11} = \frac{2}{3}\mu\gamma^2``

``\mathbf{\sigma}_{22} = -\frac{1}{3}\mu\gamma^2``

``\mathbf{\sigma}_{33} = \mathbf{\sigma_{22}}``

``\mathbf{\sigma}_{12} = \mu\gamma``

All other components are zero

### Verification
Here is a comparison of an analytic solution to the uniaxial stress boundary value problem in displacement control.
```@example
using ConstitutiveModels
using Plots

function neohookean_simple_shear()
    inputs = Dict(
        "density"         => 1.0,
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = Hyperelastic(NeoHookean())
    motion = SimpleShear(t -> t)
    out = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)

    props = initialize_props(model, inputs)
    ∇us = map(x -> x.kinematics, out)
    γs = map(x -> x[1, 2], ∇us)
    σs = map(x -> x.material_output, out)
    Zs = map(x -> x.state, out)

    μ = props[3]
    σ_11s_an = (2. / 3.) * μ * γs.^2
    σ_22s_an = -(1. / 3.) * μ * γs.^2
    σ_12s_an = μ * γs

    plot(motion, ∇us, σs, Zs, σ_11s_an, σ_22s_an, σ_12s_an)
end
neohookean_simple_shear()
```

## Uniaxial Strain
### Analytic solution
``\mathbf{\sigma}_{11} = \frac{1}{2}\kappa\left(\lambda - \frac{1}{\lambda}\right) + \frac{2}{3}\mu\left(\lambda^2 - 1\right)\lambda^{-5/3}``

``\mathbf{\sigma}_{22} = \frac{1}{2}\kappa\left(\lambda - \frac{1}{\lambda}\right) - \frac{1}{3}\mu\left(\lambda^2 - 1\right)\lambda^{-5/3}``

``\mathbf{\sigma}_{33} = \mathbf{\sigma_{22}}``

All other components are zero.

### Verification
Here is a comparison of an analytic solution to the uniaxial stress boundary value problem in displacement control.
```@example
using ConstitutiveModels
using Plots

function neohookean_uniaxial_strain()
    inputs = Dict(
        "density"         => 1.0,
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = Hyperelastic(NeoHookean())
    motion = UniaxialStrain(t -> 1 + 3t)
    out = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)

    props = initialize_props(model, inputs)
    ∇us = map(x -> x.kinematics, out)
    λs = map(x -> x[1, 1] + 1, ∇us)
    σs = map(x -> x.material_output, out)
    Zs = map(x -> x.state, out)

    κ, μ = props[2], props[3]

    σ_11s_an = 0.5 * κ * (λs .- 1 ./ λs) + 
               (2. / 3.) * μ * (λs.^2 .- 1.) .* λs .^ (-5. / 3.) 
    σ_22s_an = 0.5 * κ * (λs .- 1 ./ λs) -
               (1. / 3.) * μ * (λs.^2 .- 1.) .* λs .^ (-5. / 3.) 

    plot(motion, ∇us, σs, Zs, σ_11s_an, σ_22s_an)
end
neohookean_uniaxial_strain()
```

## Uniaxial Stress
### Analytic solution
``\mathbf{\sigma}_{11} = \frac{1}{2}\kappa\left(\lambda - \frac{1}{\lambda}\right) + \frac{2}{3}\mu\left(\lambda^2 - 1\right)\lambda^{-5/3}``

``\mathbf{\sigma}_{22} = \frac{1}{2}\kappa\left(\lambda - \frac{1}{\lambda}\right) - \frac{1}{3}\mu\left(\lambda^2 - 1\right)\lambda^{-5/3}``

``\mathbf{\sigma}_{33} = \mathbf{\sigma_{22}}``

All other components are zero.

### Verification
Here is a comparison of an analytic solution to the uniaxial stress boundary value problem in displacement control.
```@example
using ConstitutiveModels
using Plots

function neohookean_uniaxial_stress()
    inputs = Dict(
        "density"         => 1.0,
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = Hyperelastic(NeoHookean())
    motion = UniaxialStressDisplacementControl(t -> 1 + 3t)
    out = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)

    props = initialize_props(model, inputs)
    ∇us = map(x -> x.kinematics, out)
    λs = map(x -> x[1, 1] + 1, ∇us)
    σs = map(x -> x.material_output, out)
    Zs = map(x -> x.state, out)

    κ, μ = props[2], props[3]

    σ_11s_an = 0.5 * κ * (λs .- 1 ./ λs) + 
               (2. / 3.) * μ * (λs.^2 .- 1.) .* λs .^ (-5. / 3.) 
    σ_22s_an = 0.5 * κ * (λs .- 1 ./ λs) -
               (1. / 3.) * μ * (λs.^2 .- 1.) .* λs .^ (-5. / 3.) 

    plot(motion, ∇us, σs, Zs)#, σ_11s_an, σ_22s_an)
end
neohookean_uniaxial_stress()
```
