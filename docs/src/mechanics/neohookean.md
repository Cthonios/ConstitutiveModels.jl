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
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = NeoHookean()
    props = initialize_props(model, inputs)
    Δt = 0.0
    θ = 0.0
    Z = initialize_state(model)

    λs = LinRange(1., 1.25, 101)
    motion = PureShearStrain{Float64}()

    ∇us, σs, Zs = simulate_material_point(cauchy_stress, model, props, Δt, θ, Z, motion, λs)

    μ = props[2]
    σ_11s_an = (μ / 3.) * (0.5 * (λs.^2 .+ λs.^(-2)) .- 1.)
    σ_33s_an = (μ / 3.) * (2. .- λs.^2 .+ λs.^(-2))
    σ_12s_an = (μ / 2.) * (λs.^2 .- λs.^(-2))

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
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = NeoHookean()
    props = initialize_props(model, inputs)
    Δt = 0.0
    θ = 0.0
    Z = initialize_state(model)

    γs = LinRange(0., 1., 101)
    motion = SimpleShear{Float64}()

    ∇us, σs, Zs = simulate_material_point(cauchy_stress, model, props, Δt, θ, Z, motion, γs)

    μ = props[2]
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
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = NeoHookean()
    props = initialize_props(model, inputs)
    Δt = 0.0
    θ = 0.0
    Z = initialize_state(model)

    λs = LinRange(1., 4., 101)
    motion = UniaxialStrain{Float64}()

    ∇us, σs, Zs = simulate_material_point(cauchy_stress, model, props, Δt, θ, Z, motion, λs)

    κ, μ = props[1], props[2]

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

function neohookean_uniaxial_strain()
    inputs = Dict(
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = NeoHookean()
    props = initialize_props(model, inputs)
    Δt = 0.0
    θ = 0.0
    Z = initialize_state(model)

    λs = LinRange(1., 4., 101)
    motion = UniaxialStressDisplacementControl{Float64}()

    ∇us, σs, Zs = simulate_material_point(cauchy_stress, model, props, Δt, θ, Z, motion, λs)

    κ, μ = props[1], props[2]

    σ_11s_an = 0.5 * κ * (λs .- 1 ./ λs) + 
               (2. / 3.) * μ * (λs.^2 .- 1.) .* λs .^ (-5. / 3.) 
    σ_22s_an = 0.5 * κ * (λs .- 1 ./ λs) -
               (1. / 3.) * μ * (λs.^2 .- 1.) .* λs .^ (-5. / 3.) 

    plot(motion, ∇us, σs, Zs)#, σ_11s_an, σ_22s_an)
end
neohookean_uniaxial_strain()
```