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
test_finite_def_j2_tangent_large_increment()
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

function test_finite_def_j2_stress_and_tangent()
    # pk1_stress_and_material_tangent returns, from one return map, the stress
    # and the tangent of pk1_stress and material_tangent, and writes the same
    # internal variables, in the elastic range and in a plastic increment.
    model = FiniteDefJ2Plasticity()
    inputs = Dict("density" => 2700.0, "Young's modulus" => 70.0e9, "Poisson's ratio" => 0.36,
                  "yield stress" => 250.0e6, "hardening modulus" => 0.7e9)
    props = initialize_props(model, inputs)
    for ∇u in (Tensor{2, 3, Float64, 9}((1e-4, 0.0, 0.0, 0.0, -3e-5, 0.0, 0.0, 0.0, -3e-5)),
               Tensor{2, 3, Float64, 9}((0.25, 0.05, -0.2, 0.3, -0.15, 0.1, -0.1, 0.2, 0.12)))
        Z0 = initialize_state(model)
        Z1 = copy(Z0); Z2 = copy(Z0); Z3 = copy(Z0)
        P  = pk1_stress(model, props, Z0, Z1, 0.0, ∇u, 0.0)
        A  = material_tangent(model, props, Z0, Z2, 0.0, ∇u, 0.0)
        Pc, Ac = pk1_stress_and_material_tangent(model, props, Z0, Z3, 0.0, ∇u, 0.0)
        @test Pc == P
        @test Ac == A
        @test Z3 == Z1
    end
    # The default for any other model: the two functions.
    nh = Hyperelastic(NeoHookean())
    props_nh = initialize_props(nh, Dict("density" => 1.0, "Young's modulus" => 1.0,
                                         "Poisson's ratio" => 0.3))
    ∇u = Tensor{2, 3, Float64, 9}((0.1, 0.02, 0.0, -0.03, 0.05, 0.01, 0.0, 0.02, -0.04))
    Z = initialize_state(nh)
    Pc, Ac = pk1_stress_and_material_tangent(nh, props_nh, Z, copy(Z), 0.0, ∇u, 0.0)
    @test Pc == pk1_stress(nh, props_nh, Z, copy(Z), 0.0, ∇u, 0.0)
    @test Ac == material_tangent(nh, props_nh, Z, copy(Z), 0.0, ∇u, 0.0)
end

function test_finite_def_j2_tangent_large_increment()
    # BOX 9.2 returns the major-symmetric part of the Jacobian of the stress
    # update.  At a large plastic increment with hardening, the coefficient β₂
    # must use the effective shear modulus μ̄ = μ tr(b̄ᵉ_trial)/3; with μ in
    # its place the symmetric part was off by 1e-2 at an equivalent plastic
    # strain of 0.6, and exact at small strain, where μ̄ = μ.  The check: the
    # symmetric part of the central-difference Jacobian equals the tangent,
    # and what remains of the difference is the antisymmetric part.
    model = FiniteDefJ2Plasticity()
    ∇u = Tensor{2, 3, Float64, 9}((0.25, 0.05, -0.2, 0.3, -0.15, 0.1, -0.1, 0.2, 0.12))
    for H in (0.7e9, 7.0e9)
        inputs = Dict(
            "density"           => 2700.0,
            "Young's modulus"   => 70.0e9,
            "Poisson's ratio"   => 0.36,
            "yield stress"      => 250.0e6,
            "hardening modulus" => H,
        )
        props = initialize_props(model, inputs)
        Z = copy(initialize_state(model))
        pk1_stress(model, props, initialize_state(model), Z, 0.0, ∇u, 0.0)
        @test Z[10] > 0.25                      # one increment, well into the plastic range
        A = material_tangent(model, props, initialize_state(model), copy(initialize_state(model)), 0.0, ∇u, 0.0)
        h = 1.0e-6
        Afd = Tensor{4, 3}((i, j, k, l) -> begin
            δ = Tensor{2, 3}((a, b) -> (a == k && b == l) ? 1.0 : 0.0)
            Pp = pk1_stress(model, props, initialize_state(model), copy(initialize_state(model)), 0.0, ∇u + h * δ, 0.0)
            Pm = pk1_stress(model, props, initialize_state(model), copy(initialize_state(model)), 0.0, ∇u - h * δ, 0.0)
            (Pp[i, j] - Pm[i, j]) / 2h
        end)
        Am = reshape(collect(A.data), 9, 9)
        Af = reshape(collect(Afd.data), 9, 9)
        sym_err  = norm((Af + Af') / 2 - Am) / norm(Am)
        asym     = norm((Af - Af') / 2) / norm(Af)
        @test sym_err < 1.0e-7
        @test norm(Af - Am) / norm(Af) < asym + 1.0e-7
    end
end

function test_finite_def_j2_volumetric_isochoric_split()
    # The model splits as W = κ/2 (J-1)² + W_iso(b̄ᵉ) with p = κ(J-1).  The
    # isochoric functions must return W_iso, P_iso = s F⁻ᵀ and its tangent,
    # update the state exactly as the full functions do, and be independent
    # of J.  Each full quantity must equal its isochoric part plus the
    # volumetric part κ(J-1) J F⁻ᵀ (energy: κ/2 (J-1)²).
    model = FiniteDefJ2Plasticity()
    inputs = Dict(
        "density"           => 2700.0,
        "Young's modulus"   => 70.0e9,
        "Poisson's ratio"   => 0.36,
        "yield stress"      => 250.0e6,
        "hardening modulus" => 0.7e9,
    )
    props = initialize_props(model, inputs)
    λ, μ = props[2], props[3]
    @test has_volumetric_isochoric_split(model)
    @test bulk_modulus(model, props) ≈ λ + 2μ / 3
    κ = bulk_modulus(model, props)
    @test volumetric_strain(model, 1.0) == 0.0
    @test volumetric_strain_derivative(model, 1.0) == 1.0
    @test volumetric_strain_second_derivative(model, 1.3) == 0.0
    @test volumetric_strain(model, 1.3) ≈ 0.3

    # a distorted state with rotation, shear and dilatation; γ = 0.002 is
    # elastic, γ = 0.05 yields
    for γ in (2.0e-3, 5.0e-2)
        ∇u = Tensor{2, 3}((i, j) -> 0.4γ * (i + 2j) / 5 + (i == j ? 0.3γ : 0.0) + (i == 1 && j == 2 ? γ : 0.0))
        F = ∇u + one(∇u)
        J = det(F)
        Z0 = initialize_state(model)

        Za = copy(Z0); Zb = copy(Z0)
        W  = helmholtz_free_energy(model, props, Z0, Za, 0.0, ∇u, 0.0)
        Wi = isochoric_helmholtz_free_energy(model, props, Z0, Zb, 0.0, ∇u, 0.0)
        @test W ≈ Wi + κ / 2 * (J - 1)^2
        @test Za == Zb

        Za = copy(Z0); Zb = copy(Z0)
        P  = pk1_stress(model, props, Z0, Za, 0.0, ∇u, 0.0)
        Pi = isochoric_pk1_stress(model, props, Z0, Zb, 0.0, ∇u, 0.0)
        Pv = κ * (J - 1) * J * inv(F)'
        @test P ≈ Pi + Pv
        @test Za == Zb
        @test (γ > 1e-2) == (Za[10] > 0)      # the second state yields

        # the isochoric stress is deviatoric in the Kirchhoff sense: tr(P_iso Fᵀ) = 0
        @test abs(tr(dot(Pi, transpose(F)))) < 1e-8 * norm(Pi)

        # the isochoric response does not change under a pure dilatation
        c = 1.07
        Zc = copy(Z0)
        Wc = isochoric_helmholtz_free_energy(model, props, Z0, Zc, 0.0, c * F - one(F), 0.0)
        @test Wc ≈ Wi rtol = 1e-12

        # tangent: A = A_iso + ∂P_vol/∂∇u, the latter by central differences
        Za = copy(Z0); Zb = copy(Z0)
        A  = material_tangent(model, props, Z0, Za, 0.0, ∇u, 0.0)
        Ai = isochoric_material_tangent(model, props, Z0, Zb, 0.0, ∇u, 0.0)
        Pvol(g) = (Fg = g + one(g); Jg = det(Fg); κ * (Jg - 1) * Jg * inv(Fg)')
        h = 1.0e-6
        for k in 1:3, l in 1:3
            E = Tensor{2, 3}((i, j) -> (i == k && j == l) ? h : 0.0)
            D = (Pvol(∇u + E) - Pvol(∇u - E)) / (2h)
            for i in 1:3, j in 1:3
                @test isapprox(A[i, j, k, l] - Ai[i, j, k, l], D[i, j];
                               rtol = 1e-6, atol = 1e-6 * κ)
            end
        end
    end
end

test_finite_def_j2_interface()
test_finite_def_j2_uniaxial_strain()
test_finite_def_j2_cauchy_from_pk1()
test_finite_def_j2_tangent_vs_fd()
test_finite_def_j2_volumetric_isochoric_split()
test_finite_def_j2_stress_and_tangent()
