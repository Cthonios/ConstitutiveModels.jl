function test_bad_property()
    inputs = Dict("bad property" => 1.0)
    @test_throws CM.PropertyNotFoundError CM.get_property(inputs, "Young's modulus")
end

# function test_bad_

function test_default_property()
    inputs = Dict("bad property" => 1.0)
    val = CM.get_property(inputs, "Young's modulus", 70e9)
    @test val ≈ 70e9
end

function test_get_property()
    inputs = Dict(
        "Young's modulus"                  => 1.0,
        "Poisson's ratio"                  => 0.3,
        "coefficient of thermal expansion" => 0.001,
        "reference temperature"            => 60.0,
        "specific heat capacity"           => 1.0e-3,
        "thermal conductivity"             => 1.0e-4
    )
    val = CM.get_property(inputs, "reference temperature")
    @test val ≈ 60.0
end

function test_get_property_with_distribution()
    inputs = Dict(
        "Young's modulus"                  => 1.0,
        "Poisson's ratio"                  => 0.3,
        "coefficient of thermal expansion" => 0.001,
        "reference temperature"            => Dict(
            "mean"     => 60.0,
            "variance" => 1.0
        ),
        "specific heat capacity"           => 1.0e-3,
        "thermal conductivity"             => 1.0e-4
    )
    val = CM.get_property(inputs, "reference temperature")
    @test val != 60.0
end

function test_properties()
    test_bad_property()
    test_default_property()
    test_get_property()
    test_get_property_with_distribution()
end

test_properties()
