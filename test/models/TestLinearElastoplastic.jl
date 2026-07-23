function test_linear_elasto_plasticity_uniaxial_stress(model, inputs)
    props = initialize_props(model, inputs)
    λ_func = @piecewise_linear begin
        0.0, 1.00
        1.0, 1.05
    end
    motion = UniaxialStressDisplacementControl(λ_func)
    history = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)
    εs = map(x -> x.kinematics, history)
    σs = map(x -> x.material_output, history)
    Zs = map(x -> x.state, history)
    ε_xxs = map(x -> x[1, 1], εs)
    
    λ, μ = props[2], props[3]
    σ_y, H = props[4], props[5]
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
        "density"           => 1e3,
        "Young's modulus"   => E,
        "Poisson's ratio"   => ν,
        "yield stress"      => σ_y,
        "hardening modulus" => H
    )
    model = LinearElastoplastic(LinearIsotropicHardening())
    test_linear_elasto_plasticity_uniaxial_stress(model, inputs)
end

test_linear_elasto_plasticity()
