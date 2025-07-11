abstract type AbstractPlasticityModel{
    NP, 
    NS,
    E <: AbstractMechanicalModel, # elastic model
    Y <: AbstractYieldSurface
} <: AbstractMechanicalModel{NP, NS} end

function elastic_properties(model::AbstractPlasticityModel, props::SVector)
    return unpack_props(model.elastic_model, props, 1)
end

# TODO do the same as above but for yield surface
function yield_surface_properties(model::AbstractPlasticityModel, props::SVector)
    return unpack_props(model.yield_surface, props, 3)
end

function initialize_props(model::AbstractPlasticityModel, inputs::Dict{String})
    elastic_props = initialize_props(model.elastic_model, inputs)
    yield_surf_props = initialize_props(model.yield_surface, inputs)
    return vcat(elastic_props, yield_surf_props)
end

include("FeFp.jl")
include("LinearElastoPlasticity.jl")
