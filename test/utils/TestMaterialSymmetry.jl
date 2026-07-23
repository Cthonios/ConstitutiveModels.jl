function test_isotropy()
    inputs = Dict(
        "thermal conductivity" => 1e-3
    )
    symmetry = Isotropy{2}()
    props = initialize_props(symmetry, inputs, ["thermal conductivity"])
    @test all(props ≈ [1e-3])
    A = as_tensor(symmetry, props)
    @test A ≈ 1e-3 * one(SymmetricTensor{2, 3, Float64, 6})

    inputs = Dict(
        "Lamé's first constant" => 5.0,
        "shear modulus"         => 12.0
    )
    symmetry = Isotropy{4}()
    props = initialize_props(symmetry, inputs, ["Lamé's first constant", "shear modulus"])
    @test all(props ≈ [5.0, 12.0])
    C = as_tensor(symmetry, props)
    I = one(SymmetricTensor{2, 3, Float64, 6})
    @test C ≈ 5.0 * otimes(I, I) + 12.0 * (otimesu(I, I) + otimesl(I, I))
end

function test_material_symmetry()
    test_isotropy()
end

test_material_symmetry()
