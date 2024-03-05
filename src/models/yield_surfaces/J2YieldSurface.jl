struct J2YieldSurface <: YieldSurface{1, 7}
end

function J2YieldSurface(_)
  return J2YieldSurface()
end

function initialize_props(::J2YieldSurface, inputs::Dict{Symbol, Any})
  σ_y = inputs[Symbol("yield stress")]
  # return initialize_props((σ_y,))
  return SVector{1, Float64}((σ_y,))
end

function effective_stress(::J2YieldSurface, σ::M) where M <: AbstractArray{<:Number, 2}
  σ_vm = norm(dev(σ))
  return σ_vm
end
