function _test_ad_equal_analytic_for_hyper(model, props, Δt, ∇us, θ, Z_old, Z_new)
    As_ad = material_tangent.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,), (ConstitutiveModels.ForwardDiffAD(),))
    As_an = material_tangent.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(map((x, y) -> x ≈ y, As_ad, As_an))
end

#########################################################
# Gent
#########################################################

function test_gent_simple_shear(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = SimpleShear()
    γs = LinRange(0.0, 0.05, 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, γs
    )
    κ, μ, Jm = props[1], props[2], props[3]
    σ_xx_an = (2. / 3.) * Jm * μ * γs.^2 ./ (Jm .- γs.^2)
    σ_yy_an = -(1. / 3.) * Jm * μ * γs.^2 ./ (Jm .- γs.^2)
    σ_xy_an = Jm * μ * γs ./ (Jm .- γs.^2)
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_xy_an)
end

function test_gent_uniaxial_strain(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStrain()
    λs = LinRange(0.25, 4., 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    κ, μ, Jm = props[1], props[2], props[3]
    σ_xx_an = 0.5 * κ .* (λs .- 1. ./ λs) -
              (2. / 3.) * Jm * μ .* (λs.^2 .- 1.) ./
              (λs.^3 - (Jm + 3) * λs.^(5. / 3.) + 2. * λs)
    σ_yy_an = 0.5 * κ .* (λs .- 1. ./ λs) +
              (1. / 3.) * Jm * μ .* (λs.^2 .- 1.) ./
              (λs.^3 - (Jm + 3) * λs.^(5. / 3.) + 2. * λs)
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an; rtol=1e-7)
end

function test_gent()
    inputs = Dict(
        "Young's modulus" => 1.0e6,
        "Poisson's ratio" => 0.3,
        "Jm"              => 13.125
    )
    model = Gent()
    test_gent_simple_shear(model, inputs)
    test_gent_uniaxial_strain(model, inputs)
end

#########################################################
# Hencky
#########################################################
function test_hencky_uniaxial_strain(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStrain()
    λs = LinRange(0.95, 1.05, 11)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    # εs = λs .- 1.
    εs = log.(λs)
    κ, μ = props[1], props[2]
    λ = κ - (2. / 3.) * μ
    # σ_xx_an = (λ * εs .+ 2. * μ * εs) ./ λs
    # σ_xx_an = (κ * εs .+ (4. / 3.) * μ * εs) / λs
    # σ_yy_an = λ * εs
    σ_xx_an = (λ * εs .+ 2 * μ * εs) ./ λs
    σ_yy_an = (λ * εs) ./ λs
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)
end

function test_hencky()
    inputs = Dict(
        "Young's modulus" => 70.0e9,
        "Poisson's ratio" => 0.3
    )
    model = Hencky()
    # test_linear_elastic_simple_shear(model, inputs)
    # test_linear_elastic_uniaxial_strain(model, inputs)
    # test_linear_elastic_uniaxial_stress(model, inputs)
end

#########################################################
# LinearElastic
#########################################################

function test_linear_elastic_simple_shear(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = SimpleShear()
    γs = LinRange(0.0, 0.05, 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, γs
    )
    λ, μ = props[1], props[2]
    σ_xx_an = 0. * γs
    σ_yy_an = 0. * γs
    σ_xy_an = μ * γs
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_xy_an)
end

function test_linear_elastic_uniaxial_strain(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStrain()
    λs = LinRange(0.95, 1.05, 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    εs = λs .- 1.
    λ, μ = props[1], props[2]
    σ_xx_an = λ * εs .+ 2. * μ * εs
    σ_yy_an = λ * εs
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)
end

function test_linear_elastic_uniaxial_stress(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStressDisplacementControl()
    λs = LinRange(0.95, 1.05, 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    εs = map(symmetric, ∇us)
    λ, μ = props[1], props[2]
    E = μ * (3. * λ + 2. * μ) / (λ + μ)
    ν = λ / (2. * (λ + μ))
    ε_xx_an = map(x -> x[1, 1] / E, σs)
    ε_yy_an = -ν * ε_xx_an
    test_strain_eq(motion, εs, ε_xx_an, ε_yy_an)
end

function test_linear_elastic()
    inputs = Dict(
        "Young's modulus" => 70.0e9,
        "Poisson's ratio" => 0.3
    )
    model = LinearElastic()
    test_linear_elastic_simple_shear(model, inputs)
    test_linear_elastic_uniaxial_strain(model, inputs)
    test_linear_elastic_uniaxial_stress(model, inputs)
end

#########################################################
# NeoHookean
#########################################################

function test_neohookean_pure_shear_strain(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    λs = LinRange(1., 1.25, 101)
    motion = PureShearStrain{Float64}()

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    κ, μ = props[1], props[2]

    ψs_an = 0.5 .* μ .* (λs.^2 .+ λs.^(-2) .- 2)
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    σ_xx_an = (μ / 3) * (0.5 * (λs.^2 .+ λs.^(-2)) .- 1)
    σ_yy_an = (μ / 3) * (0.5 * (λs.^2 .+ λs.^(-2)) .- 1)
    σ_zz_an = (μ / 3) * (2 .- λs.^2 .- λs.^(-2))
    σ_xy_an = (μ / 2) * (λs.^2 .- λs.^(-2))
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_zz_an, σ_xy_an)

    _test_ad_equal_analytic_for_hyper(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_neohookean_simple_shear(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = SimpleShear()
    γs = LinRange(0.0, 0.5, 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, γs
    )
    κ, μ = props[1], props[2]

    ψs_an = 0.5 * μ * γs.^2
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    σ_xx_an = (2. / 3.) * μ * γs.^2
    σ_yy_an = -(1. / 3.) * μ * γs.^2
    σ_xy_an = μ * γs
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_xy_an)

    _test_ad_equal_analytic_for_hyper(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_neohookean_uniaxial_strain(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStrain()
    λs = LinRange(0.25, 4., 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    κ, μ = props[1], props[2]

    ψs_an = 0.5 * κ * (0.5 * (λs .* λs .- 1) .- log.(λs)) + 0.5 * μ * (λs.^(-2. / 3.) .* (λs.^2 .+ 2) .- 3)
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    σ_xx_an = 0.5 * κ .* (λs .- 1. ./ λs) +
              (2. / 3.) * μ .* (λs.^2 .- 1.) .* λs .^ (-5. / 3.)
    σ_yy_an = 0.5 * κ .* (λs .- 1. ./ λs) -
              (1. / 3.) * μ .* (λs.^2 .- 1.) .* λs .^ (-5. / 3.)
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)

    _test_ad_equal_analytic_for_hyper(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_neohookean()
    inputs = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3
    )
    model = NeoHookean()
    test_neohookean_pure_shear_strain(model, inputs)
    test_neohookean_simple_shear(model, inputs)
    test_neohookean_uniaxial_strain(model, inputs)
end

#########################################################
# Saint Venant-Kirchoff
#########################################################

function test_saint_venant_kirchoff_simple_shear(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = SimpleShear()
    γs = LinRange(0.0, 0.5, 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, γs
    )
    λ, μ = props[1], props[2]

    ψs_an = 0.5 .* μ .* γs.^2 .+
            0.25 .* μ .* γs.^4 .+
            0.125 .* λ .* γs.^4
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    σ_xx_an = (2μ + 0.5λ) .* γs.^2 .+ (μ + 0.5λ) .* γs.^4
    σ_yy_an = (μ + 0.5λ) .* γs.^2
    σ_zz_an = 0.5λ .* γs.^2
    σ_xy_an = μ .* γs .+ 0.5 .* (λ .+ 2μ) .* γs.^3
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_zz_an, σ_xy_an)

    _test_ad_equal_analytic_for_hyper(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_saint_venant_kirchoff_uniaxial_strain(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStrain()
    λs = LinRange(0.25, 4., 101)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )
    λ, μ = props[1], props[2]

    ψs_an = 0.125 .* λ .* (λs.^2 .- 1).^2 .+
            0.25  .* μ .* (λs.^2 .- 1).^2
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    σ_xx_an = 0.5 .* (λ + 2μ) .* λs .* (λs.^2 .- 1)
    σ_yy_an = 0.5 .* λ .* (λs.^2 .- 1) ./ λs
    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)

    _test_ad_equal_analytic_for_hyper(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_saint_venant_kirchoff()
    inputs = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3
    )
    model = SaintVenantKirchhoff()
    test_saint_venant_kirchoff_simple_shear(model, inputs)
    test_saint_venant_kirchoff_uniaxial_strain(model, inputs)
end

function test_hyperelastic_models()
    test_gent()
    test_hencky()
    test_linear_elastic()
    test_neohookean()
    test_saint_venant_kirchoff()
end

test_hyperelastic_models()
