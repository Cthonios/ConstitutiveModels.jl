# FiniteDefJ2Plasticity had no tests, which is how it came to be silently
# dropped from the module: the modularity refactor removed its `include` and
# `export` and nothing failed.  These tests pin the interface it must satisfy
# (so a future refactor that strands it again fails here) as well as its
# physics.

function test_finite_def_j2_interface()
    model = FiniteDefJ2Plasticity()

    # Reachable through the module, not just as a file on disk.
    @test isdefined(ConstitutiveModels, :FiniteDefJ2Plasticity)
    @test model isa CM.AbstractConstitutiveModel

    inputs = Dict(
        "density"           => 2700.0,
        "Young's modulus"   => 70.0e9,
        "Poisson's ratio"   => 0.36,
        "yield stress"      => 250.0e6,
        "hardening modulus" => 0.7e9,
    )
    props = initialize_props(model, inputs)

    # Density is props[1] for every model; the mass matrix and wave-speed
    # estimates of downstream codes depend on it.
    @test length(props) == num_properties(model)
    @test num_properties(model) == 5
    @test props[1] == 2700.0
    @test density(model, props, nothing, nothing, 0.0, nothing, 0.0) == 2700.0

    # p-wave modulus is λ + 2μ, read from the density-offset property layout.
    @test p_wave_modulus(model, props) ≈ props[2] + 2 * props[3]

    # State: Fᵖ = I₃ then α.
    Z = initialize_state(model)
    @test length(Z) == num_state_variables(model) == 10
    @test Z[1:9] == [1, 0, 0, 0, 1, 0, 0, 0, 1]
    @test Z[10] == 0
    @test last(state_variable_names(model)) == "eqps"
end

function test_finite_def_j2_uniaxial_strain()
    model = FiniteDefJ2Plasticity()
    inputs = Dict(
        "density"           => 2700.0,
        "Young's modulus"   => 70.0e9,
        "Poisson's ratio"   => 0.36,
        "yield stress"      => 250.0e6,
        "hardening modulus" => 0.7e9,
    )
    props = initialize_props(model, inputs)
    λ, μ, σ_y = props[2], props[3], props[4]

    # Uniaxial strain: ∇u = ε e₁⊗e₁.  Below yield the response is linear with
    # the constrained modulus λ + 2μ and no plastic flow.
    ε = 1.0e-3
    ∇u = Tensor{2, 3}((i, j) -> (i == 1 && j == 1) ? ε : 0.0)
    Z_old = initialize_state(model)
    Z_new = copy(Z_old)
    σ = cauchy_stress(model, props, Z_old, Z_new, 0.0, ∇u, 0.0)

    # von Mises check: ‖s‖ = √(3/2)·(4/3)με is still below √(2/3)σ_y here, so
    # this strain must be elastic.
    s_xx = (4 / 3) * μ * ε
    @test sqrt(1.5) * s_xx < sqrt(2 / 3) * σ_y
    @test Z_new[10] == 0.0                       # no plastic flow
    @test σ[1, 1] ≈ (λ + 2 * μ) * ε rtol = 1e-2

    # Push well past yield: plastic strain must accumulate monotonically.
    eqps = Float64[]
    for ε in (5.0e-3, 2.0e-2, 5.0e-2)
        ∇u = Tensor{2, 3}((i, j) -> (i == 1 && j == 1) ? ε : 0.0)
        Z_new = copy(initialize_state(model))
        pk1_stress(model, props, initialize_state(model), Z_new, 0.0, ∇u, 0.0)
        push!(eqps, Z_new[10])
    end
    @test all(eqps .> 0.0)
    @test issorted(eqps)
end

function test_finite_def_j2_cauchy_from_pk1()
    # σ = J⁻¹ P Fᵀ.  FiniteDefJ2Plasticity subtypes AbstractConstitutiveModel
    # directly, so it does not inherit the hyperelastic cauchy_stress fallback
    # and needs its own -- assert the two stay consistent.
    model = FiniteDefJ2Plasticity()
    inputs = Dict(
        "density"           => 2700.0,
        "Young's modulus"   => 70.0e9,
        "Poisson's ratio"   => 0.36,
        "yield stress"      => 250.0e6,
        "hardening modulus" => 0.7e9,
    )
    props = initialize_props(model, inputs)

    for ε in (1.0e-3, 3.0e-2)
        ∇u = Tensor{2, 3}((i, j) -> (i == 1 && j == 1) ? ε : 0.0)
        F  = ∇u + one(∇u)
        Zp = copy(initialize_state(model))
        Zs = copy(initialize_state(model))
        P  = pk1_stress(model, props, initialize_state(model), Zp, 0.0, ∇u, 0.0)
        σ  = cauchy_stress(model, props, initialize_state(model), Zs, 0.0, ∇u, 0.0)
        @test σ ≈ (1 / det(F)) * dot(P, transpose(F))
    end
end

function test_finite_def_j2_tangent_vs_fd()
    # The analytic tangent must be the Jacobian of the same pk1_stress the
    # residual evaluates.  This also pins the argument ORDER: the method was
    # stranded with `(props, Δt, Z_old, Z_new, ...)` while the interface is
    # `(props, Z_old, Z_new, Δt, ...)`, which would silently pass Δt as Z_old.
    model = FiniteDefJ2Plasticity()
    inputs = Dict(
        "density"           => 2700.0,
        "Young's modulus"   => 70.0e9,
        "Poisson's ratio"   => 0.36,
        "yield stress"      => 250.0e6,
        "hardening modulus" => 0.7e9,
    )
    props = initialize_props(model, inputs)

    for ε in (1.0e-3, 2.0e-2, 5.0e-2)     # elastic, then two plastic states
        ∇u = Tensor{2, 3}((i, j) -> (i == 1 && j == 1) ? ε : 0.0)

        Zs = copy(initialize_state(model))
        A  = material_tangent(model, props, initialize_state(model), Zs, 0.0, ∇u, 0.0)

        h  = 1.0e-8
        Z0 = copy(initialize_state(model))
        P0 = pk1_stress(model, props, initialize_state(model), Z0, 0.0, ∇u, 0.0)
        for k in 1:3, l in 1:3
            hh = max(h * abs(∇u[k, l]), h)
            ∇p = Tensor{2, 3}((i, j) -> ∇u[i, j] + (i == k && j == l ? hh : 0.0))
            Zp = copy(initialize_state(model))
            Pp = pk1_stress(model, props, initialize_state(model), Zp, 0.0, ∇p, 0.0)
            for i in 1:3, j in 1:3
                @test isapprox(A[i, j, k, l], (Pp[i, j] - P0[i, j]) / hh;
                               rtol = 1e-5, atol = 1e-3 * maximum(abs, P0))
            end
        end
    end
end

test_finite_def_j2_interface()
test_finite_def_j2_uniaxial_strain()
test_finite_def_j2_cauchy_from_pk1()
test_finite_def_j2_tangent_vs_fd()
