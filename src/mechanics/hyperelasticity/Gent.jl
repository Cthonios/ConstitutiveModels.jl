"""
$(TYPEDEF)
"""
struct Gent <: AbstractHyperelasticModel{3, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::Gent, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    Jm = inputs["Jm"]
    return Properties{3, eltype(elastic_props)}(
        elastic_props.κ, elastic_props.μ, Jm
    )
end

"""
``\\psi = \\frac{1}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
        - \\frac{1}{2}\\mu J_m\\ln\\left(1 - \\frac{\\bar{I}_1 - 3}{Jm}\\right)``
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::Gent,
    props, Δt,
    ∇u, θ, Z
)
    # unpack properties
    κ, μ, Jm = props[1], props[2], props[3]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
    J       = det(F)
    I_1_bar = tr(NaNMath.pow(J, -2. / 3.) * tdot(F))

    # constitutive
    ψ_vol = 0.5 * κ * (0.5 * (J^2 - 1) - NaNMath.log(J))
    ψ_dev = -μ * Jm / 2. * NaNMath.log(1. - (I_1_bar - 3.) / Jm)
    ψ     = ψ_vol + ψ_dev
    Z     = typeof(Z)()
    return ψ, Z
end