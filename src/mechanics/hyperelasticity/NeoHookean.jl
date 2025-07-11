"""
$(TYPEDEF)
"""
struct NeoHookean <: AbstractHyperelasticModel{2, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::NeoHookean, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    return Properties{2, eltype(elastic_props)}(
        elastic_props.κ, elastic_props.μ
    )
end

"""
``\\psi = \\frac{1}{2}\\left[\\frac{1}{2}\\left(J^2 - 1\\right) - \\ln J\\right]
        + \\frac{1}{2}\\mu\\left(\\bar{I}_1 - 3\\right)``
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::NeoHookean,
    props, Δt,
    ∇u, θ, Z
)
    # unpack properties
    κ, μ = props[1], props[2]

    # kinematics
    I       = one(typeof(∇u))
    F       = ∇u + I
    J       = det(F)
    I_1_bar = tr(NaNMath.pow(J, -2. / 3.) * tdot(F))

    # constitutive
    W_vol = 0.5 * κ * (0.5 * (J^2 - 1) - NaNMath.log(J))
    W_dev = 0.5 * μ * (I_1_bar - 3.)
    ψ     = W_vol + W_dev
    Z     = typeof(Z)()
    return ψ, Z
end

function pk1_stress(
    ::NeoHookean, 
    props, Δt, 
    ∇u, θ, Z
  )
  
    κ, μ    = props[1], props[2]
    F       = ∇u + one(typeof(∇u))
    J       = det(F)
    # J_23    = J^(-2. / 3.)
    J_23    = NaNMath.pow(J, -2. / 3.)
    I_1     = tr(tdot(F))
    F_inv_T = inv(F)'
    P       = 0.5 * κ * (J^2 - 1.) * F_inv_T + 
            μ * J_23 * (F - (1. / 3.) * I_1 * F_inv_T)
  
    # dummy state
    Z = typeof(Z)()
  
    return P, Z
  end