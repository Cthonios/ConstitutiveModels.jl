@muladd begin

function cos_of_acos_divided_by_3(x)
    
  x2 = x*x;
  x4 = x2*x2;

  numer = 0.866025403784438713 + 2.12714890259493060 * x + 
          ((1.89202064815951569  + 0.739603278343401613 * x) * x2 + 
            (0.121973926953064794 + x * (0.00655637626263929360 + 0.0000390884982780803443 * x)) *x4)

  denom = 1.0 + 2.26376989330935617* x + 
          ((1.80461009751278976 + 0.603976798217196003 * x) * x2 +
            (0.0783255761115461708 + 0.00268525944538021629 * x) * x4)
  
  return numer / denom
end

# takes in two methods
# the first is the gradient of the method
# used in _matrix_function
# the second is a relative method that can
# be customized to specific methods, e.g. sqrt, log, etc.
@inline function _relative_diff(f1, f2, λ1, λ2)
  λ2_safe = ifelse(λ2 == λ1, λ2 + 1.0, λ2)
  return ifelse(λ2 == λ1, f1(λ1), f2(λ1, λ2_safe))
end

@inline function eigen_sym33_non_unit(tensor::SymmetricTensor{2, 3, T, 6}) where T <: Number
  cxx = tensor[1, 1]
  cyy = tensor[2, 2]
  czz = tensor[3, 3]
  cxy = tensor[1, 2]
  cyz = tensor[2, 3]
  czx = tensor[1, 3]

  c1 = (cxx + cyy + czz) / 3.
  cxx = cxx - c1
  cyy = cyy - c1
  czz = czz - c1

  cxy_cxy = cxy * cxy
  cyz_cyz = cyz * cyz
  czx_czx = czx * czx
  cxx_cyy = cxx * cyy

  c2 = cxx_cyy + cyy * czz + czz * cxx - cxy_cxy - cyz_cyz - czx_czx
  
  c2Negative = c2 < zero(T)
  denom = ifelse(c2Negative, c2, one(T))
  ThreeOverA = ifelse(c2Negative, -3.0 / denom, one(T))
  sqrtThreeOverA = ifelse(c2Negative, sqrt(ThreeOverA), one(T))

  c3 = cxx * cyz_cyz + cyy * czx_czx - 2.0 * cxy * cyz * czx + czz * (cxy_cxy - cxx_cyy)

  rr = -0.5 * c3 * ThreeOverA * sqrtThreeOverA

  arg = min(abs(rr), one(T))

  cos_thd3 = cos_of_acos_divided_by_3(arg)

  two_cos_thd3 = 2.0 * cos_thd3 * sign(rr)

  eval2 = ifelse(c2Negative, two_cos_thd3 / sqrtThreeOverA, one(T))

  crow0 = Vec{3, T}((cxx - eval2, cxy, czx))
  crow1 = Vec{3, T}((cxy, cyy - eval2, cyz))
  crow2 = Vec{3, T}((czx, cyz, czz - eval2))

  #
  # do QR decomposition with column pivoting
  #
  k0 = crow0[1] * crow0[1] + cxy_cxy             + czx_czx
  k1 = cxy_cxy             + crow1[2] * crow1[2] + cyz_cyz
  k2 = czx_czx             + cyz_cyz             + crow2[3] * crow2[3]

  # returns zero or nan
  k0gk1 = k1 <= k0
  k0gk2 = k2 <= k0
  k1gk2 = k2 <= k1
    
  k0_largest = k0gk1 && k0gk2
  k1_largest = k1gk2 && !k0gk1
  k2_largest = !(k0_largest || k1_largest)

  k_row1_0 = ifelse(k0_largest, crow0[1], zero(T)) +
             ifelse(k1_largest, crow1[1], zero(T)) + 
             ifelse(k2_largest, crow2[1], zero(T))

  k_row1_1 = ifelse(k0_largest, crow0[2], zero(T)) +
             ifelse(k1_largest, crow1[2], zero(T)) + 
             ifelse(k2_largest, crow2[2], zero(T))

  k_row1_2 = ifelse(k0_largest, crow0[3], zero(T)) +
             ifelse(k1_largest, crow1[3], zero(T)) + 
             ifelse(k2_largest, crow2[3], zero(T))

  k_row1 = Vec{3, T}((k_row1_0, k_row1_1, k_row1_2))

  row2_0 = ifelse(k0_largest, crow1[1], crow0[1])
  row2_1 = ifelse(k0_largest, crow1[2], crow0[2])
  row2_2 = ifelse(k0_largest, crow1[3], crow0[3])

  row2 = Vec{3, T}((row2_0, row2_1, row2_2))

  row3_0 = ifelse(k2_largest, crow1[1], crow2[1])
  row3_1 = ifelse(k2_largest, crow1[2], crow2[2])
  row3_2 = ifelse(k2_largest, crow1[3], crow2[3])

  row3 = Vec{3, T}((row3_0, row3_1, row3_2))

  ki_ki = one(T) / (
    ifelse(k0_largest, k0, zero(T)) +
    ifelse(k1_largest, k1, zero(T)) +
    ifelse(k2_largest, k2, zero(T))
  )

  ki_dpr1 = ki_ki * (k_row1[1] * row2[1] + k_row1[2] * row2[2] + k_row1[3] * row2[3])
  ki_dpr2 = ki_ki * (k_row1[1] * row3[1] + k_row1[2] * row3[2] + k_row1[3] * row3[3])

  row2 = row2 - ki_dpr1 * k_row1
  row3 = row3 - ki_dpr2 * k_row1

  a0 = row2[1] * row2[1] + row2[2] * row2[2] + row2[3] * row2[3]
  a1 = row3[1] * row3[1] + row3[2] * row3[2] + row3[3] * row3[3]

  a0lea1 = a0 <= a1

  a_row2 = ifelse(a0lea1, row3, row2)
  ai_ai = one(T) / ifelse(a0lea1, a1, a0)

  evec2 = Vec{3, T}((
    k_row1[2] * a_row2[3] - k_row1[3] * a_row2[2],
    k_row1[3] * a_row2[1] - k_row1[1] * a_row2[3],
    k_row1[1] * a_row2[2] - k_row1[2] * a_row2[1]
  ))

  k_atr11 = cxx * k_row1[1] + cxy * k_row1[2] + czx * k_row1[3]
  k_atr21 = cxy * k_row1[1] + cyy * k_row1[2] + cyz * k_row1[3]
  k_atr31 = czx * k_row1[1] + cyz * k_row1[2] + czz * k_row1[3]

  a_atr12 = cxx * a_row2[1] + cxy * a_row2[2] + czx * a_row2[3]
  a_atr22 = cxy * a_row2[1] + cyy * a_row2[2] + cyz * a_row2[3]
  a_atr32 = czx * a_row2[1] + cyz * a_row2[2] + czz * a_row2[3]

  rm2xx     = (k_row1[1] * k_atr11 + k_row1[2] * k_atr21 + k_row1[3] * k_atr31) * ki_ki
  k_a_rm2xy = (k_row1[1] * a_atr12 + k_row1[2] * a_atr22 + k_row1[3] * a_atr32)
  rm2yy     = (a_row2[1] * a_atr12 + a_row2[2] * a_atr22 + a_row2[3] * a_atr32) * ai_ai
  rm2xy_rm2xy = k_a_rm2xy * k_a_rm2xy * ai_ai * ki_ki

  #
  # Wilkinson shift
  #
  b = 0.5 * (rm2xx - rm2yy)

  sqrtTerm = sqrt(b * b + rm2xy_rm2xy) * sign(b)

  eval0 = rm2yy + b - sqrtTerm
  eval1 = rm2xx + rm2yy - eval0

  rm2xx = rm2xx - eval0
  rm2yy = rm2yy - eval0

  rm2xx2 = rm2xx * rm2xx
  rm2yy2 = rm2yy * rm2yy

  fac1 = ifelse(rm2xx2 < rm2yy2, k_a_rm2xy * ai_ai, rm2xx)
  fac2 = ifelse(rm2xx2 < rm2yy2, rm2yy, ki_ki * k_a_rm2xy)

  evec0 = fac1 * a_row2 - fac2 * k_row1

  rm2xx2iszero = rm2xx2 ≈ zero(T)
  rm2xy_rm2xyiszero = rm2xy_rm2xy ≈ zero(T)
  both_zero = rm2xx2iszero & rm2xy_rm2xyiszero

  # check degeneracy
  evec0 = ifelse(both_zero, a_row2, evec0)

  evec1 = Vec{3, T}((
    evec2[2] * evec0[3] - evec2[3] * evec0[2],
    evec2[3] * evec0[1] - evec2[1] * evec0[3],
    evec2[1] * evec0[2] - evec2[2] * evec0[1]
  ))

  eval0 = eval0 + c1
  eval1 = eval1 + c1
  eval2 = eval2 + c1

  c2tol = (c1 * c1) * (-1.e-30)
  c2lsmall_neg = c2 < c2tol

  eval0 = ifelse(c2lsmall_neg, eval0, c1)
  eval1 = ifelse(c2lsmall_neg, eval1, c1)
  eval2 = ifelse(c2lsmall_neg, eval2, c1)

  evec0 = ifelse(c2lsmall_neg, evec0, basevec(Vec{3, T}, 1))
  evec1 = ifelse(c2lsmall_neg, evec1, basevec(Vec{3, T}, 2))
  evec2 = ifelse(c2lsmall_neg, evec2, basevec(Vec{3, T}, 3))

  evals = Vec{3, T}((eval0, eval1, eval2))
  evecs = Tensor{2, 3, T, 9}((
    evec0[1], evec0[2], evec0[3],
    evec1[1], evec1[2], evec1[3],
    evec2[1], evec2[2], evec2[3]
  ))

  return evals, evecs
end

@inline function eigen_sym33_unit(tensor::SymmetricTensor{2, 3, T, 6}) where T <: Number
  cmax = norm(tensor, Inf)
  cmaxInv = ifelse(cmax > zero(T), one(T) / cmax, one(T))
  scaledTensor = cmaxInv * tensor

  evals, evecs = eigen_sym33_non_unit(scaledTensor)
  evals = cmax * evals

  evec0 = evecs[:, 1]
  evec1 = evecs[:, 2]
  evec2 = evecs[:, 3]

  # evals, evec0, evec1, evec2 = eigen_sym33_non_unit(scaledTensor)
  # evals = cmax * evals

  evec0 = evec0 / norm(evec0)
  evec1 = evec1 / norm(evec1)
  evec2 = evec2 / norm(evec2)

  evecs = Tensor{2, 3, T, 9}((
    evec0[1], evec0[2], evec0[3],
    evec1[1], evec1[2], evec1[3],
    evec2[1], evec2[2], evec2[3]
  ))

  return evals, evecs
end

@inline function _matrix_function(f, A::SymmetricTensor{2, 3, T, 6}) where T <: Number
  evals, evecs = eigen_sym33_unit(A)
  new_evals = diagm(SymmetricTensor{2, 3, T, 6}, Tensors._map(f, evals))
  return dot(evecs, dot(new_evals, evecs')) |> symmetric
end

# @inline function _dmatrix_function(df, rel_diff, A::SymmetricTensor{2, 3, T, 6}) where T <: Number
@inline function _dmatrix_function(f, df, A::SymmetricTensor{2, 3, T, 6}) where T <: Number  
  evals, evecs = eigen_sym33_unit(A)
  new_evals = diagm(SymmetricTensor{2, 3, T, 6}, Tensors._map(f, evals))
  new_devals = diagm(SymmetricTensor{2, 3, T, 6}, Tensors._map(df, evals))
  fA = dot(evecs, dot(new_evals, evecs')) |> symmetric
  dfA_temp = dot(evecs, dot(new_devals, evecs')) |> symmetric
  I = one(dfA_temp)
  dfA = 0.5 * (otimesu(dfA_temp, I) + otimesl(dfA_temp, I)) |> symmetric
  return fA, dfA
end

@inline _dlog(x) = one(x) / x
# @inline _dpow(x, n) = n * NaNMath.pow(x, n - 1)
@inline _dsqrt(x) = NaNMath.sqrt(x) / 2

@inline function log_safe(A)
  if A == one(A)
    return zero(A)
  else
    return _matrix_function(NaNMath.log, A)
  end
end

# for some reason, typing this makes the derivative not register
@inline function sqrt_safe(A)
  return _matrix_function(NaNMath.sqrt, A)
end

# Define known derivative
@inline function dlog_safe(A::SymmetricTensor{2, D, T, N}) where{D, T <: Number, N}
  logA, dlogA = _dmatrix_function(NaNMath.log, _dlog, A)
  return logA, dlogA
end

@inline function dsqrt_safe(A::SymmetricTensor{2, D, T, N}) where{D, T <: Number, N}
  sqrtA, dsqrtA = _dmatrix_function(NaNMath.sqrt, _dsqrt, A)
  return sqrtA, dsqrtA
end

# Implement known derivative
@implement_gradient log_safe dlog_safe
@implement_gradient sqrt_safe dsqrt_safe

end # muladd
