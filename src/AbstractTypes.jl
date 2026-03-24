# Additional things for AD
abstract type AbstractADType end
struct ForwardDiffAD <: AbstractADType 
end
struct EnzymeAD <: AbstractADType
end

# Differentiation between
# Lagrangian (material) vs. Eulerian (spatial)
# TODO what should we call them moving forward
abstract type AbstractFrame end
struct EulerianFrame
end
struct LagrangianFrame
end

abstract type AbstractConstitutiveBase{NP, NS} end
num_properties(::AbstractConstitutiveBase{NP, NS}) where {NP, NS} = NP
num_state_variables(::AbstractConstitutiveBase{NP, NS}) where {NP, NS} = NS

# default method
function initialize_model(
    type::Type{<:AbstractConstitutiveBase}, ::Dict{String}
)
    return type()
end

"""
All models are assumed to take in a type
tag for the model and a Dict of props

This method we likely can't statically compile
due to the Any
"""
function initialize_props(
    model, 
    inputs::Dict{String}
)
    new_inputs = Dict{Symbol, Any}()
    for (key, val) in inputs
      new_inputs[Symbol(key)] = val
    end
    return initialize_props(model, new_inputs)
end

"""
Default state constructor to just return zeros
"""
function initialize_state(
    ::AbstractConstitutiveBase{NP, NS},
    float_type = Float64
) where {NP, NS}
    return zeros(float_type, NS)
end

@inline function unpack_props(
    ::AbstractConstitutiveBase{NP, NS},
    props,
    start_index::Int
# ) where {NP1, NP2, NS, T <: Number}
) where {NP, NS}
    # indices = SVector{NP1, typeof(NP1)}(start_index:start_index + NP1 - 1)
    indices = start_index:(start_index + NP - 1)
    return SVector{NP, eltype(props)}(@views props[indices])
end

abstract type AbstractConstitutiveModule{NP, NS} <: AbstractConstitutiveBase{NP, NS} end
abstract type AbstractMechanicalModule{NP, NS} <: AbstractConstitutiveModule{NP, NS} end
abstract type AbstractThermalModule{NP, NS} <: AbstractConstitutiveModule{NP, NS} end

abstract type AbstractConduction{NP, NS} <: AbstractThermalModule{NP, NS} end
# TODO add abstract interface for conduction module here
function heat_flux(::AbstractConduction, props, ∇u, θ, ∇θ)
    @assert "Implement me!"
end
"""
$(TYPEDEF)
All methods have the following type signature
(::AbstractHyperelastic, props, F::DeformationGradient, θ)

Models that call this module are responsible for first calculating
the deformation gradient from the displacement gradient
"""
abstract type AbstractHyperelastic{NP, NS} <: AbstractMechanicalModule{NP, NS} end
# needs to be implemented, no fallback
function strain_energy_density end
# defaults
function cauchy_stress(model::AbstractHyperelastic, props, F, θ)
    J = det(F)
    P = pk1_stress(model, props, F, θ)
    return (1 / J) * symmetric(P ⋅ transpose(F))
end
function material_tangent(model::AbstractHyperelastic, props, F, θ)
    return material_tangent(model, props, F, θ, ForwardDiffAD())
end
function material_tangent(model::AbstractHyperelastic, props, F, θ, ::ForwardDiffAD)
    return Tensors.gradient(z -> pk1_stress(model, props, z, θ), F)
end
function pk1_stress(model::AbstractHyperelastic, props, F, θ)
    return pk1_stress(model, props, F, θ, ForwardDiffAD())
end
function pk1_stress(model::AbstractHyperelastic, props, F, θ, ::ForwardDiffAD)
    return Tensors.gradient(z -> strain_energy_density(model, props, z, θ), F)
end
function p_wave_modulus(::AbstractHyperelastic, props)
    @assert "Implement me!"
end
abstract type AbstractIsotropicHardening{NP, NS} <: AbstractMechanicalModule{NP, NS} end
# TODO add abstract interface for isotropic hardening here
abstract type AbstractViscosity{NP, NS} <: AbstractMechanicalModule{NP, NS} end
# TODO add abstract interface for abstract viscosity here
abstract type AbstractViscoelasticity{
    NP,
    NS,
    H <: AbstractHyperelastic,
    V <: AbstractViscosity
} <: AbstractMechanicalModule{NP, NS} end
# TODO add abstract interface for abstract viscoelasticity here
abstract type AbstractYieldSurface{
    NP, 
    NS,
    IH <: AbstractIsotropicHardening
} <: AbstractConstitutiveModule{NP, NS} end
function effective_stress(::AbstractYieldSurface, ::SymmetricTensor{2, 3, T, 6}) where T <: Number
    @assert "Needs to be implemented"
end

abstract type AbstractConstitutiveModel{NP, NS} <: AbstractConstitutiveBase{NP, NS} end
# default
function density(::AbstractConstitutiveModel, props, Δt, Z_old, Z_new, args...)
    return props[1]
end
"""
Assumes it has the following structure

NameOfModelStruct:
    parameter 1: val_1
    parameter 2: val_2
    etc...
"""
function initialize(inputs::Dict{String}, model_name::String)
    parameters = inputs[model_name]
    model = initialize_model(eval(Symbol(model_name)), parameters)
    props = initialize_props(model, parameters)
    state_old = initialize_state(model)
    state_new = initialize_state(model)
    return model, props, state_old, state_new
end

# currently assumes that there is one property density
# plus the property from the mechanical modules
abstract type AbstractMechanicalModel{NP, NS} <: AbstractConstitutiveModel{NP, NS} end
function p_wave_modulus(::AbstractMechanicalModel, props)
    @assert "Implement me!"
end
abstract type AbstractThermalModel{NP, NS} <: AbstractConstitutiveModel{NP, NS} end
abstract type AbstractPlasticityModel{
    NP, 
    NS,
    E <: AbstractHyperelastic, # elastic model
    Y <: AbstractYieldSurface
} <: AbstractMechanicalModel{NP, NS} end

function elastic_properties(model::AbstractPlasticityModel, props)
    return unpack_props(model.elastic_model, props, 1)
end

# TODO do the same as above but for yield surface
function yield_surface_properties(model::AbstractPlasticityModel, props)
    return unpack_props(model.yield_surface, props, 3)
end

function initialize_props(model::AbstractPlasticityModel, inputs::Dict{String})
    elastic_props = initialize_props(model.elastic_model, inputs)
    yield_surf_props = initialize_props(model.yield_surface, inputs)
    return vcat(elastic_props, yield_surf_props)
end

abstract type AbstractThermoMechanicalModel{
    NP,
    NS,
    C <: AbstractConduction
} <: AbstractConstitutiveModel{NP, NS} end

function material_hessian(
    model::AbstractThermoMechanicalModel,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    d2ψd∇ud∇u, _ = material_tangent(model, props, Δt, ∇u, θ, Z, args...)
    d2ψdθdθ, _ = Tensors.hessian(z -> helmholtz_free_energy(
        model, props, Δt, ∇u, z, Z, args...
    ), θ)
    d2ψd∇udθ, _ = Tensors.gradient(y -> Tensors.gradient(z -> 
        helmholtz_free_energy(model, props, Δt, y, z, Z, args...), θ
    ), ∇u)

    A = tovoigt(SArray, d2ψd∇ud∇u)
    B = tovoigt(SArray, d2ψd∇udθ)

    H = vcat(
        hcat(A, B),
        hcat(B', d2ψdθdθ)
    )
    return H, Z
end


# function unpack_state(
#     ::AbstractConstitutiveModel{NP, NS1},
#     state::SVector{NS2, T},
#     start_index::Int
# ) where {NP, NS1, NS2, T <: Number}
#     indices = SVector{NS1, typeof(NS1)}(start_index:start_index + NS1 - 1)
#     return state[indices]
# end

# other then these methods...
# nothing really makes sense in general
# e.g. diffusion doesn't really have
# a "helmholtz_free_energy"
# so other interfaces will be defined
# for different physics types

# however, the following convention will be followed
# for a physics type with a method called for example
# ```physics_qoi```, we will have the following method
# signature convention
# ```physics_qoi(model, props, Δt, ∇u_1, ∇u_2, ...,  θ, Z, args...)```
# where model is the instance of a given ```AbstractConstitutiveModel```
# ```props``` is a pre-initialized flat array of properties
# ```Δt``` is the time step, ```∇u_1```, etc. 
# are the field gradients (this convention will break for hessian theories),
# ```θ``` is the absolute temperature, and ```Z``` is the
# collection of current, e.g. old, state variables.
# 
# each of these methods is assumed to return a Tensor
# e.g. rank 0, 1, 2, 4, etc., and a new set of updated
# state variables
#
# Even for models without state variables, e.g. elasticity,
# this is still necessary
