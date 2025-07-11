function test_elastic_constants()
    E = 70.e9
    ν = 0.3

    inputs = Dict(
        "Young's modulus" => E,
        "Poisson's ratio" => ν
    )
    props = ConstitutiveModels.ElasticConstants(inputs)
    @show props

    κ = E / (3. * (1. - 2. * ν))
    λ = (E * ν) / (1. + ν) / (1. - 2. * ν)
    μ = E / (2. * (1. + ν))
    @test props.κ ≈ κ
    @test props.λ ≈ λ
    @test props.ν ≈ ν
    @test props.μ ≈ μ
    @test props.E ≈ E
end

test_elastic_constants()
