# TODO
# - eventually type on stress measure rate is applied to
#   for e.g. truesdell rate on cauchy vs. kirchoff.
# - type on integration scheme forward vs. backward euler, midpoint, etc.

abstract type AbstractObjectiveStressRate end

function spatial_tangent(::AbstractObjectiveStressRate, props, Δt, σ_old, L, θ)
    λ, μ = props[1], props[2]
    δ(i, j) = i == j ? 1.0 : 0.0
    f = (i, j, k, l) -> λ * δ(i, j) * δ(k, l) + 
        μ * (δ(i, k) * δ(j, l) + δ(i, l) * δ(j, k))
    C = SymmetricTensor{4, 3, Float64, 36}(f)
    return C
end

# expected interface
function cauchy_stress end
function spatial_tangent end

struct TruesdellCauchyStressRate <: AbstractObjectiveStressRate
end

function update(model::TruesdellCauchyStressRate, props, Δt, J_old, σ_old, L, θ)
    D = symmetric(L)
    C = spatial_tangent(model, props, Δt, σ_old, L, θ)
    J = J_old + Δt * J_old * tr(L) # forward euler
    # or
    # J = J_old * exp(Δt * tr(L)) # exponential map
    σ = σ_old + Δt * (
        dcontract(C, D) + 
        dot(L, σ_old) + dot(σ_old, transpose(L)) +
        tr(L) * σ_old
    )
    return J, σ
end

struct TruesdellKirchoffStressRate <: AbstractObjectiveStressRate
end

function cauchy_stress(model::TruesdellKirchoffStressRate, props, Δt, J_old, σ_old, L, θ)
    D = symmetric(L)
    C = spatial_tangent(model, props, Δt, σ_old, L, θ)
    τ_old = J_old * σ_old
    J = J_old + Δt * J_old * tr(L) # forward euler
    τ = σ_old + Δt * (
        dcontract(C, D) + 
        dot(L, τ_old) + dot(τ_old, transpose(L))
    )
    σ = τ / J
    return J, σ
end

