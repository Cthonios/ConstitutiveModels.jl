function test_bad_property()
    inputs = Dict("bad property" => 1.0)
    @test_throws CM.PropertyNotFoundError CM.get_complex_property(inputs, "Young's modulus")
    @test_throws CM.PropertyNotFoundError CM.get_real_property(inputs, "Young's modulus")
    @test_throws CM.PropertyNotFoundError CM.get_string_property(inputs, "Young's modulus")
end

function test_default_property()
    inputs = Dict{String, Any}()
    @test CM.get_complex_property(inputs, "Shear Modulus", 1.0 + 2.0im) ≈ 1.0 + 2.0im
    @test CM.get_real_property(inputs, "Young's modulus", 70e9) ≈ 70e9
    @test CM.get_string_property(inputs, "Hardening model", "voce") == "voce"
end

function test_get_property()
    inputs = Dict(
        "Young's modulus"                  => 1.0,
        "Poisson's ratio"                  => 0.3,
        "coefficient of thermal expansion" => 0.001 + 0.1im,
        "reference temperature"            => 60.0,
        "specific heat capacity"           => 1.0e-3,
        "thermal conductivity"             => 1.0e-4,
        "hardening model"                  => "voce"
    )
    @test CM.get_complex_property(inputs, "coefficient of thermal expansion") ≈ 0.001 + 0.1im
    @test CM.get_real_property(inputs, "reference temperature") ≈ 60.0
    @test CM.get_string_property(inputs, "hardening model") == "voce"
end

function test_get_property_with_distribution()
    inputs = Dict(
        "Young's modulus"                  => 1.0,
        "Poisson's ratio"                  => 0.3,
        "coefficient of thermal expansion" => Dict(
            "mean"     => 0.001 + 0.1im,
            "variance" => 0.0001 + 0.01im
        ),
        "reference temperature"            => Dict(
            "mean"     => 60.0,
            "variance" => 1.0
        ),
        "specific heat capacity"           => 1.0e-3,
        "thermal conductivity"             => 1.0e-4
    )
    # TODO
    # @test CM.get_complex_property(inputs, "coefficient of thermal expansion") != 0.001 + 0.1im
    @test CM.get_real_property(inputs, "reference temperature") != 60.0
end

function test_properties()
    test_bad_property()
    test_default_property()
    test_get_property()
    test_get_property_with_distribution()
end

test_properties()
