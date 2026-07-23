function test_vonmises_yield_surface_linear_isotropic_hardening(inputs)
    model = VonMises(LinearIsotropicHardening())
    @test CM.has_kinematic_hardening(model) == false
    @test num_properties(model) == 2
    @test num_state_variables(model) == 1
    props = initialize_props(model, inputs)
    @test all(props .≈ [200e6, 180e6])
    state = initialize_state(model)
    @test all(state .≈ zeros(1))
end

function test_vonmises_yield_surface_no_isotropic_hardening(inputs)
    model = VonMises(NoIsotropicHardening())
    @test CM.has_kinematic_hardening(model) == false
    @test num_properties(model) == 1
    @test num_state_variables(model) == 1
    props = initialize_props(model, inputs)
    @test all(props .≈ [200e6])
    state = initialize_state(model)
    @test all(state .≈ zeros(1))
end

function test_yield_surfaces()
    inputs = Dict(
        "hardening modulus" => 180e6,
        "shear modulus"     => 30e9,
        "yield stress"      => 200e6
    )
    test_vonmises_yield_surface_linear_isotropic_hardening(inputs)
    test_vonmises_yield_surface_no_isotropic_hardening(inputs)
end

test_yield_surfaces()
