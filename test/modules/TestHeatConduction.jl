function test_isotropic_fouriers_law()
    model = FouriersLaw(Isotropy{2}())
    props = Dict(
        "thermal conductivity" => 1e-4
    )
    props = initialize_props(model, props)
    @test 1 == num_properties(model)
    @test 0 == num_state_variables(model)
    @test all([1e-4] ≈ props)
end

function test_heat_conduction()
    test_isotropic_fouriers_law()
end

test_heat_conduction()
