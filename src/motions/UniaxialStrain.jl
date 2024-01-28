@doc raw"""
Provides an analytic motion for uniaxial strain

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda & 0 & 0 \\
0       & 1 & 0 \\
0       & 0 & 1
\end{bmatrix}
```
"""
struct UniaxialStrain <: SimpleMotion
end

"""
"""
deformation_gradient(::Type{UniaxialStrain}, λ::T) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 
# deformation_gradient(::Type{UniaxialStrain}, λ::T, ::Type{<:MMatrix}) where T <: Number = 
# MMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 
# deformation_gradient(::Type{UniaxialStrain}, λ::T, ::Type{<:SMatrix}) where T <: Number = 
# SMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1., 0., 0., 0., 1.)) 
# deformation_gradient(::Type{UniaxialStrain}, λ::T, type::Type{<:Matrix}) where T <: Number = 
# T[λ 0. 0.; 0. 1. 0.; 0. 0. 1.]
