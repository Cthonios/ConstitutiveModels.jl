function test_isotropic_linear_thermal_expansion()
    model = LinearThermalExpansion(Isotropy{2}())
    props = Dict(
        "bulk modulus"                     => 1e6,
        "coefficient of thermal expansion" => 1e-3,
        "reference temperature"            => 30.0,
        "shear modulus"                    => 1e6 # just need a dummy value
    )
    props = initialize_props(model, props)
    @test all(props .≈ [30.0, -3e3])
    motion = UniaxialStrain(t -> 1 + t)
    ε = linear_strain(motion, 0.5)
    θ = 60.0
    trε = tr(ε)
    ψ = -3e3 * (θ - 30.0) * trε
    σ = -3e3 * (θ - 30.0) * one(typeof(ε))
    η = 90e3 * trε
    @test 2 == num_properties(model)
    @test 0 == num_state_variables(model)
    @test ψ ≈ helmholtz_free_energy(model, props, ε, θ)
    @test σ ≈ cauchy_stress(model, props, ε, θ)
    @test η ≈ entropy(model, props, ε, θ)
end

function test_thermal_expansion()
    test_isotropic_linear_thermal_expansion()
end

test_thermal_expansion()
