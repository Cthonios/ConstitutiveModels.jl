module ConstitutiveModels

export ADMode
export NestedADMode
export ConstitutiveModel
export num_properties
export num_state_vars
export cauchy_stress
export deformation_gradient
export helmholtz_free_energy
export pk1_stress

using DocStringExtensions
using ForwardDiff
using MuladdMacro
using StaticArrays
using Tensors

abstract type ConstitutiveModel{NP, NS} end
num_properties(::ConstitutiveModel{NP, NS}) where {NP, NS} = NP
num_state_vars(::ConstitutiveModel{NP, NS}) where {NP, NS} = NS

function initialize_props end
function initialize_state end

function initialize_props(props, ::Type{T}) where T <: Union{MVector, SVector}
  return T(props)
end

function initialize_props(props::NTuple{N, Float64}, ::Type{T}) where {N, T <: Vector}
  return collect(props)
end

function initialize_state(model, ::Type{T}) where T <: Union{MVector, SVector}
  NSV = num_state_vars(model)
  return zero(T{NSV, Float64})
end

function initialize_state(model, ::Type{T}) where T <: Vector
  return zeros(Float64, num_state_vars(model))
end

# function polar_decomposition(F::Tensor{2, 3, T, 9}) where T <: Number
#   C = tdot(F)
#   result = eigen(C)
#   U = zero(SymmetricTensor{2, 3, T, 6})
#   for n in 1:3
#     v = result.vectors[:, n]
#     U = U + sqrt(result.values[n]) * otimes(v, v)
#   end
#   R = dot(F, inv(U))
#   return R, U
# end 

# # function log_strain(F::Tensor{2, 3, T, 9}) where T <: Number
# function log_strain(F)
#   C = tdot(F)
#   result = eigen(C)
#   E = zero(SymmetricTensor{2, 3, eltype(F), 6})
#   for n in 1:3
#     v = result.vectors[:, n]
#     E = E + 0.5 * log(result.values[n]) * otimes(v, v)
#   end
#   return E
# end

struct ADMode
end

struct NestedADMode
end

function cauchy_stress end
function helmholtz_free_energy end
function pk1_stress end

"""
non-AD Cauchy stress
"""
function cauchy_stress(model::ConstitutiveModel, props, F)
  J = det(F)
  P = pk1_stress(model, props, F)
  σ = (1. / J) * dot(P, inv(F)')
  return σ
end

"""
AD for Tensors type
"""
function pk1_stress(::Type{ADMode}, model, props, F::Tensor{2, 3, T, 9}) where T <: Number
  Tensors.gradient(x -> helmholtz_free_energy(model, props, x), F)
end

"""
AD for Tensors type that uses an analytic method for the pk1 stress
"""
function pk1_tangent(::Type{ADMode}, model, props, F::Tensor{2, 3, T, 9}) where T <: Number
  Tensors.gradient(x -> pk1_stress(model, props, x), F)
end

"""
Nested AD method for Tensors type
"""
function pk1_tangent(::Type{NestedADMode}, model, props, F::Tensor{2, 3, T, 9}) where T <: Number
  Tensors.hessian(x -> helmholtz_free_energy(model, props, x), F)
end

include("Eigen.jl")
include("solvers/NewtonSolver.jl")
include("models/MechanicalModels.jl")
include("motions/Motions.jl")

# stuff for exts
function helmholtz_free_energy_sensitivies end

end # module
