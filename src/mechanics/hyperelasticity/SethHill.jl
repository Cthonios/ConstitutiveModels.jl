# Seth-Hill generalized strain family.
#
# Volumetric–deviatoric split with two integer exponents m (volumetric) and n
# (deviatoric).  Recovers neo-Hookean at m=n=1 and Hencky at m=n→0.
#
# Energy:
#   W_vol = (κ/4m²) [(Jᵐ−1)² + (J⁻ᵐ−1)²]
#   W_dev = (μ/4n²) [tr(C̄²ⁿ) + tr(C̄⁻²ⁿ) − 2 tr(C̄ⁿ) − 2 tr(C̄⁻ⁿ) + 6]
#   C̄ = J⁻²/³ C  (isochoric right Cauchy-Green tensor)
#
# PK1 stress: analytical.
# Material tangent: AD fallback via the default CommonMethods dispatch.
#
# Properties layout (NP = 4): [κ, μ, Float64(m), Float64(n)]
# m and n are stored as Float64 in props but always used as Int inside functions.

"""
$(TYPEDEF)
"""
struct SethHill <: AbstractHyperelasticModel{4, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::SethHill, inputs::Dict{String})
    ec = ElasticConstants(inputs)
    m  = Int(get(inputs, "m", 1))
    n  = Int(get(inputs, "n", 1))
    return [ec.κ, ec.μ, Float64(m), Float64(n)]
end

"""
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::SethHill,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    κ, μ   = props[1], props[2]
    m, n   = Int(props[3]), Int(props[4])

    F      = ∇u + one(typeof(∇u))
    J      = det(F)
    Jm     = J^m;  Jmm = 1.0 / Jm
    C      = tdot(F)
    Cbar   = symmetric(J^(-2/3) * C)

    Cbarn  = _matrix_function(x -> x^n,   Cbar)
    Cbarmn = _matrix_function(x -> x^(-n), Cbar)
    Cbar2n = _matrix_function(x -> x^(2n),  Cbar)
    Cbarm2n= _matrix_function(x -> x^(-2n), Cbar)

    W_vol = κ / (4m^2) * ((Jm - 1)^2 + (Jmm - 1)^2)
    W_dev = μ / (4n^2) * (tr(Cbar2n) + tr(Cbarm2n) - 2tr(Cbarn) - 2tr(Cbarmn) + 6)
    return W_vol + W_dev
end

function pk1_stress(
    ::SethHill,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    κ, μ    = props[1], props[2]
    m, n    = Int(props[3]), Int(props[4])

    F       = ∇u + one(typeof(∇u))
    J       = det(F)
    Jm      = J^m;  J2m = Jm^2;  Jmm = 1.0 / Jm;  Jm2m = 1.0 / J2m

    C       = tdot(F)
    F_inv_T = inv(F)'
    Cbar    = symmetric(J^(-2/3) * C)
    Cbarmn  = _matrix_function(x -> x^(-n), Cbar)
    Cbarn   = _matrix_function(x -> x^n,    Cbar)
    Cbar2n  = _matrix_function(x -> x^(2n), Cbar)
    Cbarm2n = _matrix_function(x -> x^(-2n), Cbar)

    trCbarn   = tr(Cbarn);   trCbarmn  = tr(Cbarmn)
    trCbar2n  = tr(Cbar2n);  trCbarm2n = tr(Cbarm2n)

    P_vol = κ / (2m) * (J2m - Jm - Jm2m + Jmm) * F_inv_T

    scalar = (1 / 3) * (-trCbar2n + trCbarn + trCbarm2n - trCbarmn)
    combo  = Tensor{2, 3, eltype(Cbar2n), 9}(Cbar2n) - Tensor{2, 3, eltype(Cbarn), 9}(Cbarn) -
             Tensor{2, 3, eltype(Cbarm2n), 9}(Cbarm2n) + Tensor{2, 3, eltype(Cbarmn), 9}(Cbarmn)
    P_dev  = μ / n * (scalar * F_inv_T + F_inv_T ⋅ combo)

    return P_vol + P_dev
end

# material_tangent falls through to the AD default in CommonMethods.jl
