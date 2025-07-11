# LinearElastic
```@autodocs
Modules = [ConstitutiveModels]
Order   = [:type, :function]
Pages   = ["LinearElastic.jl"]
```

## Simple shear
### Analytic solution
Fill me out

### Verification
Here is a comparison of an analytic solution to the uniaxial stress boundary value problem in displacement control.
```@example
using ConstitutiveModels
using Plots

function linearelastic_simple_shear()
    inputs = Dict(
        "Young's modulus" => 1.0,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = LinearElastic()
    props = initialize_props(model, inputs)
    Δt = 0.0
    θ = 0.0
    Z = initialize_state(model)

    γs = LinRange(0., 0.01,101)
    motion = SimpleShear{Float64}()

    ∇us, σs, Zs = simulate_material_point(cauchy_stress, model, props, Δt, θ, Z, motion, γs)

    μ = props[2]
    σ_11s_an = 0. * γs
    σ_22s_an = 0. * γs
    σ_12s_an = μ * γs

    plot(motion, ∇us, σs, Zs, σ_11s_an, σ_22s_an, σ_12s_an)
end
linearelastic_simple_shear()
```

## Uniaxial Strain
### Analytic solution
Fill me out

### Verification
Here is a comparison of an analytic solution to the uniaxial stress boundary value problem in displacement control.
```@example
using ConstitutiveModels
using Plots

function linearelastic_uniaxial_strain()
    inputs = Dict(
        "Young's modulus" => 100e3,#u"MPa",
        "Poisson's ratio" => 0.3
    )

    model = LinearElastic()
    props = initialize_props(model, inputs)
    Δt = 0.0
    θ = 0.0
    Z = initialize_state(model)

    λs = LinRange(1., 1.001, 101)
    motion = UniaxialStrain{Float64}()

    # hardcoded parameters for simple models
    ∇us, σs, Zs = simulate_material_point(cauchy_stress, model, props, Δt, θ, Z, motion, λs)

    λ, μ = props[1], props[2]
    σ_11s_an = λ * (λs .- 1.) .+ 2. * μ * (λs .- 1.)
    σ_22s_an = λ * (λs .- 1.)
    plot(motion, ∇us, σs, Zs, σ_11s_an, σ_22s_an)

end

linearelastic_uniaxial_strain()
```