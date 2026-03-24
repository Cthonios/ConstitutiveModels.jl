abstract type AbstractThermoMechanicalModel{
    NP, NS, 
    M <: AbstractMechanicalModule,
    T <: AbstractThermalModel 
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

include("LinearThermoElastic.jl")
