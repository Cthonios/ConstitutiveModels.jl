@doc raw"""
Provides an analytic motion for simple shear.

This is
```math
\mathbf{F} = \begin{bmatrix}
1 & \gamma & 0 \\
0 & 1      & 0 \\
0 & 0      & 1
\end{bmatrix}
```
"""
struct SimpleShear <: SimpleMotion
end

"""
"""
deformation_gradient(::Type{SimpleShear}, γ::T) where T <: Number = 
Tensor{2, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))
# deformation_gradient(::Type{SimpleShear}, γ::T, ::Type{<:MMatrix}) where T <: Number = 
# MMatrix{3, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))
# deformation_gradient(::Type{SimpleShear}, γ::T, ::Type{<:SMatrix}) where T <: Number = 
# SMatrix{3, 3, T, 9}((1., 0., 0., γ, 1., 0., 0., 0., 1.))
# deformation_gradient(::Type{SimpleShear}, γ::T, type::Type{<:Matrix}) where T <: Number = 
# T[1. γ 0.; 0. 1. 0.; 0. 0. 1.]
