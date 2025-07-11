struct LinearThermoElastic{
    NP, NS,
    M <: LinearElastic,
    T <: FouriersLaw
} <: AbstractThermoMechanicalModel{NP, NS, M, T}
    elastic_model::M
    thermal_model::T
end

function LinearThermoElastic()
    elastic_model = LinearElastic()
    thermal_model = FouriersLaw()
    return LinearThermoElastic{
        6, 0, typeof(elastic_model), typeof(thermal_model)
    }(elastic_model, thermal_model)
end

function initialize_props(::LinearThermoElastic, inputs::Dict{String})
    elastic_props = ElasticConstants(inputs)
    α = inputs["coefficient of thermal expansion"]
    θ_0 = inputs["reference temperature"]
    c = inputs["specific heat capacity"]
    k = inputs["thermal conductivity"]
    β = -3. * elastic_props.κ * α
    return Properties{6, eltype(β)}(
        elastic_props.λ, elastic_props.μ,
        β, θ_0, c, k
    )
end

function helmholtz_free_energy(
    model::LinearThermoElastic,
    props, Δt,
    ∇u, θ, Z
)
    # unpack properties
    λ, μ, β, θ_0, c, _ = props
    elastic_props = SVector{2, eltype(props)}(λ, μ)

    # kinematics
    ε = symmetric(∇u)

    # purely elastic part
    ψ_e, _ = helmholtz_free_energy(
        model.elastic_model,
        elastic_props, Δt,
        ∇u, θ, Z
    )

    # purely therml part
    ψ_t = -(c / (2. * θ_0)) * (θ - θ_0)^2
    
    # mixed term
    ψ_mixed = β * (θ - θ_0) * tr(ε)

    Z = typeof(Z)()
    return ψ_e + ψ_t + ψ_mixed, Z
end
