# Push-forward of a 4th-order Lagrangian moduli tensor ℂ to the first elasticity
# tensor A = ∂P/∂F.
#
# Formula (Norma convention, validated against FD tangent):
#
#   A_{iJkL} = S_{LJ} δ_{ik} + Σ_{m,n} F_{im} ℂ_{mJLn} F_{kn}
#
# where:
#   ℂ  = 2 ∂S/∂C  (Lagrangian moduli, includes the factor-of-2 from symmetric C)
#   S  = PK2 stress
#   F  = deformation gradient
#
# Works with both Tensor{4,3} and SymmetricTensor{4,3} inputs for ℂ and S.
# Returns Tensor{4,3} (∂P/∂F has no symmetry in general).
#
# Tensors.jl stores Tensor{4,3} in column-major order:
#   linear index = i + 3*(j-1) + 9*(k-1) + 27*(l-1)
@inline function _convect_tangent(
    ℂ,                       # Tensor{4,3} or SymmetricTensor{4,3}: Lagrangian moduli
    S,                       # Tensor{2,3} or SymmetricTensor{2,3}: PK2 stress
    F::Tensor{2, 3, T, 9},   # deformation gradient
) where {T <: Number}
    data = ntuple(Val(81)) do lin
        l, rem = divrem(lin - 1, 27)
        k, rem = divrem(rem,      9)
        j, i   = divrem(rem,      3)
        i += 1; j += 1; k += 1; l += 1
        val = S[l, j] * (i == k ? one(T) : zero(T))
        for m in 1:3, n in 1:3
            val += F[i, m] * ℂ[m, j, l, n] * F[k, n]
        end
        val
    end
    return Tensor{4, 3, T, 81}(data)
end
