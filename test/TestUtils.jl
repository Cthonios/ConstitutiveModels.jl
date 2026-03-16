# this test method uses the Tensors eigen method
function _generate_random_psd_tensor()
    A = rand(Tensor{2, 3, Float64, 9})
    return tdot(A)
end

function _matrix_function_test(f, A::SymmetricTensor{2, 3, T, 6}, args...) where T <: Number
    evals, evecs = eigen(A)
    new_evals = diagm(SymmetricTensor{2, 3, T, 6}, Tensors._map(x -> f(x, args...), evals))
    return dot(evecs, dot(new_evals, evecs')) |> symmetric
end

function _dmatrix_function(f, A, args...)
    return Tensors.gradient(z -> _matrix_function_test(x -> f(x, args...), z), A)
end

function test_exp()
    A = _generate_random_psd_tensor()
    expA_test = _matrix_function_test(exp, A)
    expA = exp(A)
    @test expA_test ≈ expA

    @test log(expA) ≈ A

    dexpA_test = _dmatrix_function(exp, A)
    dexpA = ConstitutiveModels.dexp(A)
    @test dexpA_test ≈ dexpA
end

function test_log()
    A = _generate_random_psd_tensor()
    logA_test = _matrix_function_test(log, A)
    logA = log(A)
    @test logA_test ≈ logA

    @test exp(logA) ≈ A

    dlogA_test = _dmatrix_function(log, A)
    dlogA = ConstitutiveModels.dlog(A)
    @test dlogA_test ≈ dlogA
end

function test_pow()
    A = _generate_random_psd_tensor()

    # test consistency with sqrt
    powA_test = _matrix_function_test(pow, A, 0.5)
    powA = pow(A, 0.5)
    @test powA_test ≈ powA

    @test dot(powA, powA) ≈ A

    dpow_test = _dmatrix_function(pow, A, 0.5)
    dpowA = ConstitutiveModels.dpow(A, 0.5)
    @test dpow_test ≈ dpowA

    # test random power between 0 and 1
    m = rand(1)[1]
    powA_test = _matrix_function_test(pow, A, m)
    powA = pow(A, m)
    @test powA_test ≈ powA

    @test pow(powA, 1 / m) ≈ A

    dpow_test = _dmatrix_function(pow, A, m)
    dpowA = ConstitutiveModels.dpow(A, m)
    @test dpow_test ≈ dpowA
end

function test_sqrt()
    A = _generate_random_psd_tensor()
    sqrtA_test = _matrix_function_test(sqrt, A)
    sqrtA = sqrt(A)
    @test sqrtA_test ≈ sqrtA

    @test dot(sqrtA, sqrtA) ≈ A

    dsqrtA_test = _dmatrix_function(sqrt, A)
    dsqrtA = ConstitutiveModels.dsqrt(A)
    @test dsqrtA_test ≈ dsqrtA
end

function test_utils()
    test_exp()
    test_log()
    test_pow()
    test_sqrt()
end

test_utils()
