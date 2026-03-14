# Saint Venant-Kirchhoff hyperelastic model.
#
# The finite-deformation analogue of linear elasticity: the same quadratic
# energy density, but expressed in terms of the Green-Lagrange strain
# E = ½(C − I) rather than the infinitesimal strain ε = sym(∇u).
#
#   ψ = ½ λ (tr E)² + μ E : E
#   S = λ tr(E) I + 2μ E        (PK2)
#   ℂ = λ I⊗I + 2μ I⊙I         (constant Lagrangian moduli)
#
# Note: ℂ here is 2 ∂S/∂C per the Norma convention used in _convect_tangent.
# Working through ∂S/∂C = ½(λ I⊗I + 2μ I⊙I) gives ℂ = 2∂S/∂C = λ I⊗I + 2μ I⊙I,
# but because ∂E/∂C = ½ I⊙I (symmetric, not ½ factor absorbed),
# the correct formula is:
#   2 ∂S_{AB}/∂C_{CD} = λ δ_{AB} δ_{CD} + μ(δ_{AC} δ_{BD} + δ_{AD} δ_{BC})
# which is ℂ_vol + ℂ_dev below.  Verified against FD tangent.

"""
$(TYPEDEF)
"""
struct SaintVenantKirchhoff <: AbstractHyperelasticModel{2, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::SaintVenantKirchhoff, inputs::Dict{String})
    ec = ElasticConstants(inputs)
    return [ec.λ, ec.μ]
end

"""
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::SaintVenantKirchhoff,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    λ, μ = props[1], props[2]
    F    = ∇u + one(typeof(∇u))
    C    = tdot(F)
    E    = 0.5 * (C - one(C))         # Green-Lagrange strain (SymmetricTensor{2,3})
    trE  = tr(E)
    return 0.5λ * trE^2 + μ * dcontract(E, E)
end

function pk1_stress(
    ::SaintVenantKirchhoff,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    λ, μ = props[1], props[2]
    F    = ∇u + one(typeof(∇u))
    C    = tdot(F)
    E    = 0.5 * (C - one(C))
    S    = λ * tr(E) * one(E) + 2.0μ * E   # PK2
    return F ⋅ S                             # PK1 = F S
end

# Analytical material tangent ∂P/∂F.
#
# The Lagrangian moduli ℂ = 2 ∂S/∂C are CONSTANT for Saint Venant-Kirchhoff:
#   ℂ = λ I⊗I + 2μ I⊙I
# (same index structure as linear elasticity; coefficient = 1 not 2 because the
#  factor of 2 is absorbed into the ∂E/∂C = ½ I⊙I chain rule).
# Push-forward via _convect_tangent gives ∂P/∂F.
function material_tangent(
    ::SaintVenantKirchhoff,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    λ, μ = props[1], props[2]
    F    = ∇u + one(typeof(∇u))
    C    = tdot(F)
    E    = 0.5 * (C - one(C))
    S    = λ * tr(E) * one(E) + 2.0μ * E

    I2   = one(symmetric(C))                    # SymmetricTensor{2,3} identity
    ℂ    = λ * otimes(I2, I2) +
           μ * (otimesu(I2, I2) + otimesl(I2, I2))  # = λ I⊗I + 2μ I⊙I

    return _convect_tangent(ℂ, Tensor{2, 3, eltype(S), 9}(S), F)
end
