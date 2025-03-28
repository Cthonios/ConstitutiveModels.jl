@muladd begin

@inline function deformation_gradient(∇u)
  F = ∇u + one(∇u)
  # return TwoPointTensor(F)
  return F
end

@inline function green_lagrange_strain(F::Tensor)
  C = right_cauchy_green(F)
  return (1 / 2) * (C - one(C))
end

# function hencky_strain(F::TwoPointTensor{T}) where T <: Number
#   C = right_cauchy_green(F)
#   # evals, evecs = eigen(C.val)
#   # log_λs = Tensors._map(log, evals)
#   # return MaterialTensor((1 / 2) * (
#   #   log_λs[1] * otimes(evecs[:, 1]) + 
#   #   log_λs[2] * otimes(evecs[:, 2]) + 
#   #   log_λs[3] * otimes(evecs[:, 3])
#   # ))
#   evals, evecs = eigen_sym33_unit(C.val)
#   log_λs = diagm(SymmetricTensor{2, 3, T, 6}, Tensors._map(log, evals))
#   return MaterialTensor((1 / 2) * dot(evecs, dot(log_λs, evecs')) |> symmetric)
# end


function hencky_strain(C::SymmetricTensor)
  evals, evecs = eigen_sym33_unit(C)
  log_λs = diagm(SymmetricTensor{2, 3, eltype(C), 6}, Tensors._map(log, evals))
  # return MaterialTensor((1 / 2) * dot(evecs, dot(log_λs, evecs')) |> symmetric)
  return (1 / 2) * dot(evecs, dot(log_λs, evecs')) |> symmetric
end

@inline function left_cauchy_green(F::Tensor)
  return dott(F)
end

@inline function linear_strain(∇u)
  # return SpatialTensor(symmetric(∇u))
  return symmetric(∇u)
end

@inline function right_cauchy_green(F::Tensor)
  return tdot(F)
end

end
