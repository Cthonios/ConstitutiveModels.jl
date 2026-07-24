struct MyHeatConduction <: CM.AbstractHeatConduction
end

function test_bad_heat_conduction_interface()
    model = MyHeatConduction()
    props = Float64[]
    ∇u = zero(Tensor{2, 3, Float64, 9})
    θ = 30.0
    ∇θ = Vec{3, Float64}((1.0, 2.0, 3.0))
    @test_throws CM.UnImplementedMethodError dissipation(model, props, ∇u, θ, ∇θ)
    @test_throws CM.UnImplementedMethodError heat_flux(model, props, ∇u, θ, ∇θ)
end

function test_bad_fouriers_law_symmetry()
    model = FouriersLaw(MyBadSymmetry{2}())
    props = Dict(
        "thermal conductivity" => 1e-4
    )
    @test_throws AssertionError initialize_props(model, props)
end

function test_isotropic_fouriers_law_no_deformation()
    model = FouriersLaw(Isotropy{2}())
    props = Dict(
        "thermal conductivity" => 1e-4
    )
    props = initialize_props(model, props)
    @test 1 == num_properties(model)
    @test 0 == num_state_variables(model)
    @test all([1e-4] ≈ props)

    # non-linear case
    ∇u = zero(Tensor{2, 3, Float64, 9})
    θ = 30.0
    ∇θ = Vec{3, Float64}((1.0, 2.0, 3.0))
    d = dissipation(model, props, ∇u, θ, ∇θ)
    q = heat_flux(model, props, ∇u, θ, ∇θ)
    # @show d
    # @show 
    @test d ≈ 1e-4 / 900 * dot(∇θ, ∇θ)
    @test q ≈ -1e-4 * ∇θ

end

function test_heat_conduction()
    test_bad_fouriers_law_symmetry()
    test_bad_heat_conduction_interface()
    test_isotropic_fouriers_law_no_deformation()
end

test_heat_conduction()
