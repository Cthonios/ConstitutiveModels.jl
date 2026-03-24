"""
$(TYPEDEF)
"""
struct FouriersLaw{Frame <: AbstractFrame} <: AbstractConduction{1, 0}
    frame::Frame
end

"""
$(TYPEDSIGNATURES)
defaults to the Lagrangian frame unless
the frame option is present.
"""
function intialize_model(::Type{<:FouriersLaw}, inputs::Dict{String})
    if haskey(inputs, "frame")
        frame = eval(Symbol(inputs["frame"]))
    else
        frame = LagrangianFrame()
    end
    return FouriersLaw(frame)
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
function heat_flux(::FouriersLaw{<:EulerianFrame}, props, ∇u, θ, ∇θ)

    # unpack props
    k = props[1]

    # caclulate heat flux
    q = -k * ∇θ

    return q
end

"""
$(TYPEDSIGNATURES)
Calculates the heat flux for 
a simple Fourier's law and also returns an empy
state

``\\mathbf{q} = -k\\nabla\\mathbf{C}^{-1}\\theta``
"""
function heat_flux(::FouriersLaw{<:LagrangianFrame}, props, ∇u, θ, ∇θ)

    # unpack props
    k = props[1]

    # kinematics
    F = ∇u + one(∇u)
    C = tdot(F)
    Cinv = inv(C)
  
    # caclulate heat flux
    q = -k * Cinv ⋅ ∇θ

    return q
end
