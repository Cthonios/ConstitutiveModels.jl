# For AD w.r.t. to scalars, i.e. temperature
function Tensors._extract_gradient(v::Tuple{ForwardDiff.Dual, SVector}, t::Number)
    η = Tensors._extract_gradient(v[1], t)
    return η, ForwardDiff.value.(v[2])
end

# For AD w.r.t. to second order tensors, i.e. the deformation gradient
function Tensors._extract_gradient(v::Tuple{ForwardDiff.Dual, SVector}, t::Tensor{2, dim}) where {dim}
    P = Tensors._extract_gradient(v[1], t)
    return P, ForwardDiff.value.(v[2])
end
  
# For nested AD w.r.t to second order tensors
function Tensors._extract_gradient(
    v::Tuple{Tensor{2, 3, <:ForwardDiff.Dual}, SVector}, 
    t::Tensor{2, dim}
) where {dim}
    A = Tensors._extract_gradient(v[1], t)
    return A, ForwardDiff.value.(v[2])
end
  
# # For eigen decomp
# function Tensors._extract_gradient(
#     v::Tuple{Vec{dim, <:ForwardDiff.Dual}, Tensor{2, 3, <:ForwardDiff.Dual}},
#     t::SymmetricTensor{2, dim, T}
# ) where {dim, T <: Number}
#     @show "hur"
#     v1 = Tensors._extract_gradient(v[1], t)
#     display(v1)
#     @assert false
# end

# This catches the case where AD can't trace the temperature
# in hyperelastic models so we just set the entropy to zero
# since mathematically speaking dψ/dθ = 0.0 for hyperelstic
# models
# function Tensors._extract_gradient(v::Tuple{T, SVector}, ::T) where T <: Number
#     return -0.0, v[2]
# end 

# function setup_ad end