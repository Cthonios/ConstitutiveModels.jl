abstract type AbstractHeatConduction <: AbstractConstitutiveModule end

# should we rename below to something like push-forward or pull-back?
# TODO check against Gurtin-Fried-Anand
# non-linear kinematics
function map_temperature_gradient(∇u::Tensor{2, 3, T, 9}, ∇θ) where T <: Number
    F = ∇u + one(∇u)
    C_inv = inv(tdot(F))
    return C_inv * ∇θ
end

# linear kinematics
function map_temperature_gradient(::SymmetricTensor{2, 3, T, 6}, ∇θ) where T <: Number
    return ∇θ
end

num_state_variables(::AbstractHeatConduction) = 0

# expected interface

# TODO setup expected "tangent" interface.
# function conductivity end

function dissipation(::AbstractHeatConduction, props, ∇u, θ, ∇θ)
    throw(UnImplementedMethodError("Need to implement dissipation method for heat conduction module."))
end

function heat_flux(::AbstractHeatConduction, props, ∇u, θ, ∇θ)
    throw(UnImplementedMethodError("Need to implement heat_flux method for heat conduction module."))
end

# implementations below

# TODO make use MaterialSymmetry.jl interface
struct FouriersLaw{S <: AbstractMaterialSymmetry{2}} <: AbstractHeatConduction
    symmetry::S

    function FouriersLaw(symmetry = Isotropy{2}())
        new{typeof(symmetry)}(symmetry)
    end
end

"""
$(TYPEDSIGNATURES)
"""
function initialize_props(model::FouriersLaw, inputs::Dict{String})
    if isa(model.symmetry, Isotropy{2})
        k = inputs["thermal conductivity"]
        return SVector{1, typeof(k)}(k)
    else
        @assert false "Currently unsupported material symmetry type $(model.symmetry)"
    end
end

num_properties(model::FouriersLaw) = num_properties(model.symmetry)

# TODO check against Gurtin-Fried-Anand
function dissipation(model::FouriersLaw, props, kin, θ, ∇θ)
    q = heat_flux(model, props, kin, θ, ∇θ)
    return -(1 / (θ * θ)) * dot(q, ∇θ)
end

function heat_flux(::FouriersLaw, props, kin, θ, ∇θ)
    # k = props[1]
    k = as_tensor(model.symmetry, props)
    ∇θ = map_temperature_gradient(kin, ∇θ)
    # return -k * ∇θ
    return -dot(k, ∇θ)
end
