function test_linear_elasto_plasticity_uniaxial_stress(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStressDisplacementControl()
    λs = LinRange(1.0, 1.05, 101) # be careful not to start in e.g. compression

    # TODO this below is dumb
    # the test fails on stress equality
    # with cauchy stress as the method
    ∇us, σs, Zs = simulate_material_point(
        pk1_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    εs = map(symmetric, ∇us)
    ε_xxs = map(x -> x[1, 1], εs)
    
    λ, μ = props[1], props[2]
    σ_y, H = props[3], props[4]
    E = μ * (3. * λ + 2. * μ) / (λ + μ)

    ε_p_an = Float64[]
    σ_xx_an = Float64[]
    for ε in ε_xxs
        σ = E * ε
        if σ <= σ_y
            push!(ε_p_an, 0.0)
            push!(σ_xx_an, σ)
        else
            ε_p = (E * ε - σ_y) / (E + H)
            # σ = (σ_y + H * ε) / (1 + (H / E))
            σ = σ_y + H * ε_p
            push!(ε_p_an, ε_p)
            push!(σ_xx_an, σ)
        end
    end

    for (σ, Z, σ_xx, ε_p) in zip(σs, Zs, σ_xx_an, ε_p_an)
        @test σ[1, 1] ≈ σ_xx
        @test Z[7] ≈ ε_p
    end
end

function test_linear_elasto_plasticity()
    E = 70.e9
    ν = 0.3
    σ_y = 200.e6
    H = 180.e6
    inputs = Dict(
        "Young's modulus"           => E,
        "Poisson's ratio"           => ν,
        "yield surface"             => "VonMisesYieldSurface",
        "isotropic hardening model" => "LinearIsotropicHardening",
        "yield stress"              => σ_y,
        "hardening modulus"         => H
    )
    model = LinearElastoPlasticity(
        VonMisesYieldSurface(
            LinearIsotropicHardening()
        )
    )
    test_linear_elasto_plasticity_uniaxial_stress(model, inputs)
end

test_linear_elasto_plasticity()
