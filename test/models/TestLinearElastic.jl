function test_linear_elastic_model_simple_shear(model, inputs)
    props = initialize_props(model, inputs)
    γ_func = @piecewise_linear begin
        0.0, 0.00
        1.0, 0.25
    end
    motion = SimpleShear(γ_func)
    history = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)
    # ∇us = map(x -> x.kinematics, history)
    εs = map(x -> x.kinematics, history)
    σs = map(x -> x.material_output, history)
    γs = 2. * map(x -> x[1, 2], εs)
    λ, μ = props[2], props[3]
    σ_xx_an = 0. * γs
    σ_yy_an = 0. * γs
    σ_xy_an = μ * γs
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_xy_an)
end

function test_linear_elastic_model_uniaxial_strain(model, inputs)
    props = initialize_props(model, inputs)
    λ_func = @piecewise_linear begin
        0.0, 0.25
        1.0, 4.00
    end
    motion = UniaxialStrain(λ_func)
    history = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)
    εs = map(x -> x.kinematics[1, 1], history)
    σs = map(x -> x.material_output, history)
    λ, μ = props[2], props[3]
    σ_xx_an = λ * εs .+ 2. * μ * εs
    σ_yy_an = λ * εs
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)
end

function test_linear_elastic_model_uniaxial_stress(model, inputs)
    props = initialize_props(model, inputs)
    λ_func = @piecewise_linear begin
        0.0, 0.25
        1.0, 4.00
    end
    motion = UniaxialStressDisplacementControl(λ_func)
    history = simulate_material_point(cauchy_stress, model, inputs, motion, 1.0)
    εs = map(x -> x.kinematics, history)
    σs = map(x -> x.material_output, history)
    λ, μ = props[2], props[3]
    E = μ * (3. * λ + 2. * μ) / (λ + μ)
    ν = λ / (2. * (λ + μ))
    ε_xx_an = map(x -> x[1, 1] / E, σs)
    ε_yy_an = -ν * ε_xx_an
    test_strain_eq(motion, εs, ε_xx_an, ε_yy_an)
end

function test_linear_elastic_model()
    inputs = Dict(
        "density"         => 1e3,
        "Young's modulus" => 70.0e9,
        "Poisson's ratio" => 0.3
    )
    model = LinearElastic()
    test_linear_elastic_model_simple_shear(model, inputs)
    test_linear_elastic_model_uniaxial_strain(model, inputs)
    test_linear_elastic_model_uniaxial_stress(model, inputs)
end

test_linear_elastic_model()
