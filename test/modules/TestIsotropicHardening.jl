function test_linear_isotropic_hardening(inputs)
    model = LinearIsotropicHardening()
    @test num_properties(model) == 2
    @test num_state_variables(model) == 1
    props = initialize_props(model, inputs)
    @test all(props .≈ [200e6, 180e6])
    state = initialize_state(model)
    @test all(state .≈ [0.0])

    eqps = 0.3
    @test CM.hardening_energy(model, props, eqps) ≈ 200e6 * 0.3 + 90e6 * 0.09
    @test CM.hardening_gradient(model, props, eqps) ≈ 200e6 + 180e6 * 0.3
    @test CM.hardening_hessian(model, props, eqps) ≈ 180e6
    @test CM.yield_surface_radius(model, props, eqps) ≈ sqrt(2. / 3.) * (200e6 + 180e6 * 0.3)
end

function test_no_isotropic_hardening(inputs)
    model = NoIsotropicHardening()
    @test num_properties(model) == 1
    @test num_state_variables(model) == 1
    props = initialize_props(model, inputs)
    @test all(props .≈ [200e6])
    state = initialize_state(model)
    @test all(state .≈ Float64[])

    eqps = 0.3
    @test CM.hardening_energy(model, props, eqps) ≈ 200e6 * 0.3
    @test CM.hardening_gradient(model, props, eqps) ≈ 200e6
    @test CM.hardening_hessian(model, props, eqps) ≈ 0.0
    @test CM.yield_surface_radius(model, props, eqps) ≈ sqrt(2. / 3.) * 200e6
end

function test_isotropic_hardening()
    inputs = Dict(
        "hardening modulus" => 180e6,
        "yield stress"      => 200e6
    )
    test_linear_isotropic_hardening(inputs)
    test_no_isotropic_hardening(inputs)
end

test_isotropic_hardening()
