struct VonMises <: YieldSurface{1, 0}
end

function initialize_properties(::VonMises, inputs::Dict{Symbol, Any})
  σ_y = inputs[Symbol("yield stress")]
  # return initialize_props((σ_y,))
  return SVector{1, Float64}((σ_y,))
end

function effective_stress(::VonMises, σ)
  σ_vm = norm(dev(σ))
  # σ_vm = sqrt(dcontract(dev(σ), dev(σ)))
  return σ_vm
end
