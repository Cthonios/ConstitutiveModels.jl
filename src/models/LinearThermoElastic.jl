struct LinearThermoElastic{
    NP, NS,
    C <: FouriersLaw,
    M <: LinearElastic
} <: AbstractThermoMechanicalModel{NP, NS, C}
    thermal_model::C
    elastic_model::M
end

function LinearThermoElastic()
    thermal_model = FouriersLaw()
    elastic_model = LinearElastic()
    return LinearThermoElastic{
        6, 0, typeof(thermal_model), typeof(elastic_model)
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
    props, Δt, Z_old, Z_new,
    ∇u, θ
)
    # unpack properties
    λ, μ, β, θ_0, c, _ = props
    elastic_props = SVector{2, eltype(props)}(λ, μ)

    # kinematics
    ε = symmetric(∇u)

    # purely elastic part
    ψ_e, _ = helmholtz_free_energy(
        model.elastic_model,
        elastic_props, Δt, Z_old, Z_new,
        ∇u, θ
    )

    # purely therml part
    ψ_t = -(c / (2. * θ_0)) * (θ - θ_0)^2
    
    # mixed term
    ψ_mixed = β * (θ - θ_0) * tr(ε)

    return ψ_e + ψ_t + ψ_mixed
end
