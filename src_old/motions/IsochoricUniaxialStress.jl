@doc raw"""
Provides an analytic motion for uniaxial stress
assuming perfect incompressibility.

This is
```math
\mathbf{F} = \begin{bmatrix}
\lambda & 0 &                        0 \\
0       & \frac{1}{\sqrt{\lambda}} & 0 \\
0       & 0                        & \frac{1}{\sqrt{\lambda}}
\end{bmatrix}
```
"""
struct IsochoricUniaxialStress <: SimpleMotion
end

"""
"""
deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T) where T <: Number = 
Tensor{2, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))
# deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T, type::Type{<:MMatrix}) where T <: Number = 
# MMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))
# deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T, type::Type{<:SMatrix}) where T <: Number = 
# SMatrix{3, 3, T, 9}((λ, 0., 0., 0., 1. / sqrt(λ), 0., 0., 0., 1. / sqrt(λ)))
# deformation_gradient(::Type{IsochoricUniaxialStress}, λ::T, type::Type{<:Matrix}) where T <: Number = 
# T[λ 0. 0.; 0. 1. / sqrt(λ) 0.; 0. 0. 1. / sqrt(λ)]
