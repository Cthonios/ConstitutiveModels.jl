"""
$(TYPEDEF)
"""
struct FouriersLaw <: AbstractThermalModel{1, 0}
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(::FouriersLaw, inputs::Dict{String})
    k = inputs["thermal conductivity"]
    return SVector{1, typeof(k)}(k)
end

"""
$(TYPEDSIGNATURES)
Calculates the heat flux for 
a simple Fourier's law and also returns an empy
state

``\\mathbf{q} = -k\\nabla\\theta``
"""
function heat_flux(::FouriersLaw, props, Δt, ∇θ, θ, Z)

    # unpack props
    k = props[1]
  
    # caclulate heat flux
    q = -k * ∇θ
  
    # dummy state
    Z = typeof(Z)()
  
    return q, Z
end
