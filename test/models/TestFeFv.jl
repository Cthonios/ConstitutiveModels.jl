function log_symm(A)
    return 0.5 * log(tdot(A))
end

function fefv_test_inputs()
    return Dict(
        "density"                   => 1100.0,
        "bulk modulus"  => 855.0,
        "shear modulus" => 0.855,
        "G neq"                     => 5.0,
        "relaxation time"           => 25.0,
    )
end

function test_fefv_interface()
    model = FeFv()

    @test isdefined(ConstitutiveModels, :FeFv)
    @test model isa CM.AbstractHyperelasticModel

    props = initialize_props(model, fefv_test_inputs())

    # Density is props[1]; G_neq, τ hardcoded at props[4], props[5] in the
    # source implies exactly 2 Hencky properties sit between them -- pin
    # that layout so a future Hencky change doesn't silently misalign it.
    @test length(props) == num_properties(model)
    @test num_properties(model) == 5
    @test props[1] == 1100.0
    @test props[4] == 5.0
    @test props[5] == 25.0

    # Regression guard for the model.model / model.model_eq field-name and
    # flattening bugs.
    @test length(property_names(model)) == num_properties(model)

    Z = initialize_state(model)
    @test length(Z) == num_state_variables(model) == 9
    @test Z == vec(collect(one(Tensor{2, 3, Float64, 9})))

    @test length(state_variable_names(model)) == num_state_variables(model)

    @test p_wave_modulus(model, props) == p_wave_modulus(model.model_eq, props[2:3])
end

function test_fefv_uniaxial_loading()
    model = FeFv()
    props = initialize_props(model, fefv_test_inputs())
    τ = props[5]

    strain_rate = 1.0e-2
    total_time  = 100.0
    n_steps     = 100
    dt          = total_time / n_steps
    times       = dt .* (0:(n_steps - 1))

    Z_old = initialize_state(model)
    Ev_11 = zeros(n_steps)
    Ee_11 = zeros(n_steps)
    Ee_22 = zeros(n_steps)

    for (n, t) in enumerate(times)
        F = Tensor{2, 3}((i, j) -> begin
            if i == j == 1
                exp(strain_rate * t)
            elseif i == j
                1.0
            else
                0.0
            end
        end)
        ∇u = F - one(F)

        Z_new = copy(Z_old)
        pk1_stress(model, props, Z_old, Z_new, dt, ∇u, 0.0)

        Fv = ConstitutiveModels.unpack_state(model, Z_new)
        Ev = log_symm(Fv)
        Fe = dot(F, inv(Fv))
        Ee = log_symm(Fe)

        Ev_11[n] = Ev[1, 1]
        Ee_11[n] = Ee[1, 1]
        Ee_22[n] = Ee[2, 2]

        Z_old = Z_new
    end

    # Closed-form solution of the linear Maxwell-branch ODE for this
    # loading history (independent derivation from the model code).
    e_v_11_analytic = @. (2 / 3) * strain_rate * times -
                          (2 / 3) * strain_rate * τ * (1 - exp(-times / τ))
    e_e_11_analytic = @. strain_rate * times - e_v_11_analytic
    e_e_22_analytic = 0.5 .* e_v_11_analytic

    @test maximum(abs.(Ev_11 .- e_v_11_analytic)) < 1.5e-3
    @test maximum(abs.(Ee_11 .- e_e_11_analytic)) < 1.5e-3
    @test maximum(abs.(Ee_22 .- e_e_22_analytic)) < 1.5e-3
end

function test_fefv_cauchy_from_pk1()
    model = FeFv()
    props = initialize_props(model, fefv_test_inputs())
    dt = 1.0

    for ε in (1.0e-3, 3.0e-2)
        ∇u = Tensor{2, 3}((i, j) -> (i == 1 && j == 1) ? ε : 0.0)
        F = ∇u + one(∇u)

        Zp = copy(initialize_state(model))
        Zs = copy(initialize_state(model))
        P = pk1_stress(model, props, initialize_state(model), Zp, dt, ∇u, 0.0)
        σ = cauchy_stress(model, props, initialize_state(model), Zs, dt, ∇u, 0.0)

        @test σ ≈ (1 / det(F)) * dot(P, transpose(F))
    end
end

function test_fefv_tangent_vs_fd()
    model = FeFv()
    props = initialize_props(model, fefv_test_inputs())
    dt = 1.0

    # Exercise the tangent both from a fresh Fv = I state and from a state
    # with already-accumulated viscous flow: an index-order mistake in the
    # Fv_old^{-T} pull-back would be invisible at Fv_old = I (since
    # inv(I)^T = I hides the transpose/ordering) but should show up once
    # Fv_old is anisotropic.
    Z_fresh = initialize_state(model)
    Z_flowed = copy(Z_fresh)
    let ∇u_pre = Tensor{2, 3}((i, j) -> (i == 1 && j == 1) ? 0.05 : 0.0)
        pk1_stress(model, props, initialize_state(model), Z_flowed, dt, ∇u_pre, 0.0)
    end

    for Z_old in (Z_fresh, Z_flowed), ε in (1.0e-3, 2.0e-2)
        ∇u = Tensor{2, 3}((i, j) -> begin
            if i == 1 && j == 1
                ε
            elseif i == 2 && j == 1
                0.3 * ε
            else
                0.0
            end
        end)

        Zs = copy(Z_old)
        A = material_tangent(model, props, Z_old, Zs, dt, ∇u, 0.0)

        h = 1.0e-8
        Z0 = copy(Z_old)
        P0 = pk1_stress(model, props, Z_old, Z0, dt, ∇u, 0.0)
        for k in 1:3, l in 1:3
            hh = max(h * abs(∇u[k, l]), h)
            ∇p = Tensor{2, 3}((i, j) -> ∇u[i, j] + (i == k && j == l ? hh : 0.0))
            Zp = copy(Z_old)
            Pp = pk1_stress(model, props, Z_old, Zp, dt, ∇p, 0.0)
            for i in 1:3, j in 1:3
                @test isapprox(A[i, j, k, l], (Pp[i, j] - P0[i, j]) / hh;
                               rtol = 1.0e-4, atol = 1.0e-4 * maximum(abs, P0) + 1.0e-8)
            end
        end
    end
end

test_fefv_interface()
test_fefv_uniaxial_loading()
test_fefv_cauchy_from_pk1()
test_fefv_tangent_vs_fd()
