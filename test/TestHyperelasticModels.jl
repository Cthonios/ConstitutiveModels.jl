function _test_ad_equal_analytic_for_hyper_material_tangent(
    model, props, Δt, ∇us, θ, Z_old, Z_new;
    atol = 1e-10, rtol = 1e-10
)
    As_ad = material_tangent.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,), (ConstitutiveModels.ForwardDiffAD(),))
    As_an = material_tangent.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    # display(map((x, y) -> x - y, As_ad, As_an))
    @test all(map((x, y) -> isapprox(x, y, atol = atol, rtol = rtol), As_ad, As_an))
end

function _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
    Ps_ad = pk1_stress.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,), (ConstitutiveModels.ForwardDiffAD(),))
    Ps_an = pk1_stress.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    # display(map((x, y) -> x - y, Ps_ad, Ps_an))
    # @assert false
    @test all(map((x, y) -> isapprox(x, y, atol = 1e-10, rtol = 1e-10), Ps_ad, Ps_an))
end

#########################################################
# ArrudaBoyce
#########################################################
function test_arruda_boyce_uniaxial_strain(model, inputs)
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
    κ, C1, C2 = props[1], props[2], props[3]

    # # analytical Helmholtz free energy
    # J = λs
    # J_m_13 = 1.0 ./ cbrt.(J)
    # J_m_23 = J_m_13 .* J_m_13
    # J_m_43 = J_m_23 .* J_m_23
    # I1 = λs.^2 .+ 2
    # I2 = 0.5 * (I1.^2 .- λs.^4 .- 2.)
    # I1_bar = J_m_23 .* I1
    # I2_bar = J_m_23 .* J_m_23 .* 0.5 .* (I1.^2 .- (λs.^4 .+ 2))
    # ψs_an = 0.5 * κ * (0.5 * (J.^2 .- 1) .- log.(J)) .+ C1 .* (I1_bar .- 3) .+ C2 .* (I2_bar .- 3)
    
    # ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    # @test all(ψs_an .≈ ψs)

    # # TODO
    # # analytical Cauchy stress components (compressible uniaxial)
    # σ_xx_an = 0.5 * κ .* (λs .- 1. ./ λs) +
    #       (2/3) * 2C1 .* (λs.^2 .- 1.) .* λs.^(-5/3) +
    #       (4/3) * C2 .* (λs.^2 .- 1.) .* λs.^(-7/3)

    # σ_yy_an = 0.5 * κ .* (λs .- 1. ./ λs) -
    #       (1/3) * 2C1 .* (λs.^2 .- 1.) .* λs.^(-5/3) -
    #       (2/3) * C2 .* (λs.^2 .- 1.) .* λs.^(-7/3)

    # test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)

    # check AD vs analytic tangent and PK1 stress
    # _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_arruda_boyce()
    inputs = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3,
        "n"               => 25.0
    )
    model = ArrudaBoyce()
    # test_mooney_rivlin_simple_shear(model, inputs)
    # test_arruda_boyce_uniaxial_strain(model, inputs)
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

    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
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

    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
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
function test_hencky_simple_shear(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = SimpleShear()
    # AD failing on zero deformation
    γs = LinRange(0.0005, 0.05, 11)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, γs
    )
    κ, μ = props[1], props[2]
    σs_an = map(∇us) do ∇u
        F   = ∇u + one(∇u)
        J   = det(F)
        C   = tdot(F)
        IC  = inv(C)
        E   = 0.5 * log(C)
        A   = κ * tr(E) * one(C) + 2μ * dev(E)
        S   = IC ⋅ A
        σ   = symmetric((1 / J) * F ⋅ S ⋅ transpose(F))
        return σ
    end
    σ_xx_an = map(x -> x[1, 1], σs_an)
    σ_yy_an = map(x -> x[2, 2], σs_an)
    σ_zz_an = map(x -> x[3, 3], σs_an)
    σ_xy_an = map(x -> x[1, 2], σs_an)

    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_zz_an, σ_xy_an)

    _test_ad_equal_analytic_for_hyper_material_tangent(
        model, props, Δt, ∇us, θ, Z_old, Z_new;
        atol = 1e-9, rtol=1e-9
    )
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_hencky_uniaxial_strain(model, inputs)
    props = initialize_props(model, inputs)
    Z_old = initialize_state(model)
    Z_new = initialize_state(model)
    Δt = 0.0
    θ = 0.0

    motion = UniaxialStrain()
    λs = LinRange(1.0005, 1.05, 11)

    ∇us, σs, _ = simulate_material_point(
        cauchy_stress, model, props, Δt, θ, Z_old, Z_new, motion, λs
    )

    κ, μ = props[1], props[2]
    σs_an = map(∇us) do ∇u
        F   = ∇u + one(∇u)
        J   = det(F)
        C   = tdot(F)
        IC  = inv(C)
        E   = 0.5 * log(C)
        A   = κ * tr(E) * one(C) + 2μ * dev(E)
        S   = IC ⋅ A
        σ   = symmetric((1 / J) * F ⋅ S ⋅ transpose(F))
        return σ
    end
    σ_xx_an = map(x -> x[1, 1], σs_an)
    σ_yy_an = map(x -> x[2, 2], σs_an)

    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)

    _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_hencky()
    inputs = Dict(
        "Young's modulus" => 70.0e3,
        "Poisson's ratio" => 0.3
    )
    model = Hencky()
    test_hencky_simple_shear(model, inputs)
    # test_hencky_uniaxial_strain(model, inputs)
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
function test_mooney_rivlin_simple_shear(model, inputs)
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
    _, C1, C2 = props

    # analytical Helmholtz free energy (isochoric invariants: I1_bar = I2_bar = 3 + γ^2)
    ψs_an = (C1 + C2) .* γs.^2
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    # analytical Cauchy stress components
    σ_xy_an = 2 * (C1 + C2) .* γs
    σ_xx_an = (2/3) * (C1 + C2) .* γs.^2
    σ_yy_an = -(1/3) * (C1 + C2) .* γs.^2

    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an, σ_xy_an)

    _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_mooney_rivlin_uniaxial_strain(model, inputs)
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
    κ, C1, C2 = props[1], props[2], props[3]

    # analytical Helmholtz free energy
    J = λs
    J_m_13 = 1.0 ./ cbrt.(J)
    J_m_23 = J_m_13 .* J_m_13
    J_m_43 = J_m_23 .* J_m_23
    I1 = λs.^2 .+ 2
    I2 = 0.5 * (I1.^2 .- λs.^4 .- 2.)
    I1_bar = J_m_23 .* I1
    I2_bar = J_m_23 .* J_m_23 .* 0.5 .* (I1.^2 .- (λs.^4 .+ 2))
    ψs_an = 0.5 * κ * (0.5 * (J.^2 .- 1) .- log.(J)) .+ C1 .* (I1_bar .- 3) .+ C2 .* (I2_bar .- 3)
    
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    # TODO
    # analytical Cauchy stress components (compressible uniaxial)
    σ_xx_an = 0.5 * κ .* (λs .- 1. ./ λs) +
          (2/3) * 2C1 .* (λs.^2 .- 1.) .* λs.^(-5/3) +
          (4/3) * C2 .* (λs.^2 .- 1.) .* λs.^(-7/3)

    σ_yy_an = 0.5 * κ .* (λs .- 1. ./ λs) -
          (1/3) * 2C1 .* (λs.^2 .- 1.) .* λs.^(-5/3) -
          (2/3) * C2 .* (λs.^2 .- 1.) .* λs.^(-7/3)

    test_stress_eq(motion, σs, σ_xx_an, σ_yy_an)

    # check AD vs analytic tangent and PK1 stress
    # _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_mooney_rivlin()
    inputs = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3,
        "C1"              => 0.5,
        "C2"              => 0.5
    )
    model = MooneyRivlin()
    # test_mooney_rivlin_simple_shear(model, inputs)
    test_mooney_rivlin_uniaxial_strain(model, inputs)
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

    _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
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

    _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
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

    _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
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

    _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
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

    _test_ad_equal_analytic_for_hyper_material_tangent(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
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

function test_seth_hill_simple_shear(model, inputs)
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

    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_seth_hill_uniaxial_strain(model, inputs)
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
    # λ, μ = props[1], props[2]
    κ, μ, m, n = props[1], props[2], props[3], props[4]

    Js = λs
    # volumetric energy
    ψs_vol = κ / (4 * m^2) * ((Js.^m .- 1).^2 .+ (Js.^(-m) .- 1).^2)

    # isochoric right Cauchy-Green
    C_xx = λs.^2
    C_yy = 1 .^2
    C_zz = 1 .^2
    Jm23 = Js.^(-2/3)
    Cbar_xx = C_xx .* Jm23
    Cbar_yy = C_yy .* Jm23
    Cbar_zz = C_zz .* Jm23

    # Deviatoric invariants
    trCbar_n    = Cbar_xx.^n .+ Cbar_yy.^n .+ Cbar_zz.^n
    trCbar_invn = Cbar_xx.^(-n) .+ Cbar_yy.^(-n) .+ Cbar_zz.^(-n)
    trCbar_2n   = Cbar_xx.^(2n) .+ Cbar_yy.^(2n) .+ Cbar_zz.^(2n)
    trCbar_inv2 = Cbar_xx.^(-2n) .+ Cbar_yy.^(-2n) .+ Cbar_zz.^(-2n)

    ψs_dev = μ / (4 * n^2) .* (trCbar_2n .+ trCbar_inv2 .- 2*trCbar_n .- 2*trCbar_invn .+ 6)

    ψs_an = ψs_vol + ψs_dev
    ψs = helmholtz_free_energy.((model,), (props,), (Δt,), ∇us, (θ,), (Z_old,), (Z_new,))
    @test all(ψs_an .≈ ψs)

    # _test_ad_equal_analytic_for_hyper(model, props, Δt, ∇us, θ, Z_old, Z_new)
    _test_ad_equal_analytic_for_hyper_pk1_stress(model, props, Δt, ∇us, θ, Z_old, Z_new)
end

function test_seth_hill()
    inputs_neohookean = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3,
        "m"               => 1,
        "n"               => 1
    )
    inputs_1_2 = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3,
        "m"               => 1,
        "n"               => 2
    )
    inputs_2_1 = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3,
        "m"               => 2,
        "n"               => 1
    )
    inputs_2_2 = Dict(
        "Young's modulus" => 1.0,
        "Poisson's ratio" => 0.3,
        "m"               => 2,
        "n"               => 2
    )
    model = SethHill()
    test_seth_hill_simple_shear(model, inputs_neohookean)
    test_seth_hill_uniaxial_strain(model, inputs_neohookean)

    test_seth_hill_simple_shear(model, inputs_1_2)
    test_seth_hill_uniaxial_strain(model, inputs_1_2)

    test_seth_hill_simple_shear(model, inputs_2_1)
    test_seth_hill_uniaxial_strain(model, inputs_2_1)

    test_seth_hill_simple_shear(model, inputs_2_2)
    test_seth_hill_uniaxial_strain(model, inputs_2_2)
end

function test_hyperelastic_models()
    test_arruda_boyce()
    test_gent()
    test_hencky()
    test_linear_elastic()
    test_mooney_rivlin()
    test_neohookean()
    test_saint_venant_kirchoff()
    test_seth_hill()
end

test_hyperelastic_models()
