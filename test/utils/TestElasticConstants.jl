function compare(props_1, props_2)
    @test props_1.κ ≈ props_2.κ
    @test props_1.λ ≈ props_2.λ
    @test props_1.ν ≈ props_2.ν
    @test props_1.μ ≈ props_2.μ
    @test props_1.E ≈ props_2.E
end

function test_elastic_constants()
    # young's modulus/poisson's ratio combo
    E = 70.e9
    ν = 0.3

    inputs = Dict(
        "Young's modulus" => E,
        "Poisson's ratio" => ν
    )
    props = CM.ElasticConstants(inputs)
    @show props

    κ = E / (3. * (1. - 2. * ν))
    λ = (E * ν) / (1. + ν) / (1. - 2. * ν)
    μ = E / (2. * (1. + ν))
    @test props.κ ≈ κ
    @test props.λ ≈ λ
    @test props.ν ≈ ν
    @test props.μ ≈ μ
    @test props.E ≈ E

    # reuse above values for rest of tests
    # young's modulus/bulk modulus combo
    inputs = Dict(
        "Young's modulus" => 70e9,
        "bulk modulus"    => 5.833333333333333e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # young's modulus/lame constant combo
    inputs = Dict(
        "Young's modulus"       => 70e9,
        "Lamé's first constant" => 4.038461538461538e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # young's modulus/shear modulus combo
    inputs = Dict(
        "Young's modulus" => 70e9,
        "shear modulus"   => 2.6923076923076923e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # young's modulus/bad property combo
    inputs = Dict(
        "Young's modulus" => 70e9,
        "some other prop" => 1.0
    )
    @test_throws CM.ElasticConstantsError CM.ElasticConstants(inputs)

    # poisson's ratio/bulk modulus constant
    inputs = Dict(
        "Poisson's ratio" => 0.3,
        "bulk modulus"    => 5.833333333333333e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # poisson's ratio/lame constant combo
    inputs = Dict(
        "Poisson's ratio"       => 0.3,
        "Lamé's first constant" => 4.038461538461538e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # poisson's ratio/shear modulus
    inputs = Dict(
        "Poisson's ratio" => 0.3,
        "shear modulus"   => 2.6923076923076923e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # poisson's ratio/bad property combo
    inputs = Dict(
        "Poisson's ratio"  => 0.3,
        "some other props" => 1.0
    )

    @test_throws CM.ElasticConstantsError CM.ElasticConstants(inputs)

    # bulk modulus/lame param combo
    inputs = Dict(
        "bulk modulus"          => 5.833333333333333e10,
        "Lamé's first constant" => 4.038461538461538e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # bulk modulus/shear modulus
    inputs = Dict(
        "bulk modulus"  => 5.833333333333333e10,
        "shear modulus" => 2.6923076923076923e10
    )
    compare(props, CM.ElasticConstants(inputs))

    # bulk modulus/bad property combo
    inputs = Dict(
        "bulk modulus"     => 0.3,
        "some other props" => 1.0
    )
    @test_throws CM.ElasticConstantsError CM.ElasticConstants(inputs)
    
    # lame param/shear modulus combo
    inputs = Dict(
        "Lamé's first constant" => 4.038461538461538e10,
        "shear modulus"         => 2.6923076923076923e10
    )
    compare(props, CM.ElasticConstants(inputs))
end

test_elastic_constants()
