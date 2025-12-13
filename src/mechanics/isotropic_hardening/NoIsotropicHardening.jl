struct NoIsotropicHardening <: AbstractIsotropicHardening{1, 0}
end

function initialize_props(::NoIsotropicHardening, inputs::Dict{String})
    σ_y = inputs["yield stress"]
    return [σ_y]
end

energy(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray = props[1] * eqps
radius(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray = sqrt(2. / 3.) * props[1]
slope(::NoIsotropicHardening, props::V, eqps::T) where {V <: AbstractArray, T <: Number}  = zero(T)

function hardening_increment(
    model::NoIsotropicHardening,
    props, μ, σ_eff::T, α_old
) where T <: Number
  
    f = σ_eff - radius(model, props, α_old)
    
    if f <= zero(T)
      Δγ = 0.0
    else
      Δγ = f / (2. * μ)
    end
    return Δγ
end
