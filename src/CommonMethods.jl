const ModelsWithMechanics = Union{
    <:AbstractMechanicalModel,
    <:AbstractThermoMechanicalModel
}
const ModelsWithThermal = Union{
    <:AbstractThermalModel,
    <:AbstractThermoMechanicalModel
}
const ModelsWithThermoMechanical = Union{
    ModelsWithMechanics,
    ModelsWithThermal
}

"""
$(TYPEDSIGNATURES)
"""
function cauchy_stress(
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    F = ∇u + one(∇u)
    J = det(F)
    P, Z = pk1_stress(model, props, Δt, ∇u, θ, Z, args...)
    σ = (1. / J) * dot(P, F')
    return σ, Z
end

function entropy(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    return entropy(
        model, props, Δt,
        ∇u, θ, Z,
        ForwardDiffAD(),
        args...
    )
end

function heat_capacity(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    return heat_capacity(
        model, props, Δt,
        ∇u, θ, Z,
        ForwardDiffAD(),
        args...
    )
end

# function algorithmic_energy end
"""
Expected abstract interface where \n
Inputs:\n
```props```      - array of properties\n
```Δt````        - time step\n
```∇θ```         - absolute temperature gradient\n
```θ```          - absolute temperature\n
```Z```          - current (old) state variable array\n
```args```       - potentially additional arguments for special models\n
Ouputs:\n
``ψ``            - helmholtz free energy\n
``\\mathcal{Z}`` - new state variable array (different from input ```Z```)
$(TYPEDSIGNATURES)
"""
function helmholtz_free_energy(
    ::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    @assert false "Needs to be implemented"
end

"""
Expected abstract interface where \n
Inputs:\n
```props```      - array of properties\n
```Δt````        - time step\n
```∇θ```         - absolute temperature gradient\n
```θ```          - absolute temperature\n
```Z```          - current (old) state variable array\n
```args```       - potentially additional arguments for special models\n
Ouputs:\n
``\\mathbf{P}``  - first Piola-Kirchoff stress tensor\n
``\\mathcal{Z}`` - new state variable array (different from input ```Z```)

If this method is not defined for a model, it will fallback to 
the current default AD interface in ```ConstitutiveModels.jl```
by differentiating the output of ```helmholtz_free_energy``` with respect
to ∇u, e.g.

``\\mathbb{A} = \\frac{\\partial^2\\psi}{\\partial\\nabla\\mathbf{u}\\partial\\nabla\\mathbf{u}} = 
\\frac{\\partial^2\\psi}{\\partial\\mathbf{F}\\partial\\mathbf{F}}``
$(TYPEDSIGNATURES)
"""
function material_tangent(
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    return material_tangent(
        model, props, Δt,
        ∇u, θ, Z,
        ForwardDiffAD(),
        args...
    )
end

"""
Expected abstract interface where \n
Inputs:\n
```props```      - array of properties\n
```Δt````        - time step\n
```∇θ```         - absolute temperature gradient\n
```θ```          - absolute temperature\n
```Z```          - current (old) state variable array\n
```args```       - potentially additional arguments for special models\n
Ouputs:\n
``\\mathbf{P}``  - first Piola-Kirchoff stress tensor\n
``\\mathcal{Z}`` - new state variable array (different from input ```Z```)

If this method is not defined for a model, it will fallback to 
the current default AD interface in ```ConstitutiveModels.jl```
by differentiating the output of ```helmholtz_free_energy``` with respect
to ∇u, e.g.

``\\mathbf{P} = \\frac{\\partial\\psi}{\\partial\\nabla\\mathbf{u}} = 
\\frac{\\partial\\psi}{\\partial\\mathbf{F}}``
$(TYPEDSIGNATURES)
"""
function pk1_stress(
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z,
    args...
)
    return pk1_stress(
        model, props, Δt,
        ∇u, θ, Z,
        ForwardDiffAD(),
        args...
    )
end

# ForwardDiff wrappers

function entropy(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z::SVector{NS, T},
    ::ForwardDiffAD,
    args...
)::Tuple{T, SVector{NS, T}} where {NS, T}
    results = Tensors.gradient(z -> 
        helmholtz_free_energy(model, props, Δt, ∇u, z, Z, args...),
        θ
    )
    η = -results[1]
    return η, results[2]
end

function heat_capacity(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z::SVector{NS, T},
    ::ForwardDiffAD,
    args...
)::Tuple{T, SVector{NS, T}} where {NS, T}
    results = Tensors.hessian(z -> 
        helmholtz_free_energy(model, props, Δt, ∇u, z, Z, args...),
        θ
    )
    c = -θ * results[1]
    return c, results[2]
end

function material_tangent(
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z,
    ::ForwardDiffAD,
    args...
)
    return Tensors.gradient(z -> 
        pk1_stress(model, props, Δt, z, θ, Z, ForwardDiffAD(), args...),
        ∇u
    )
end

function pk1_stress(
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z,
    ::ForwardDiffAD,
    args...
)
    return Tensors.gradient(z -> 
        helmholtz_free_energy(model, props, Δt, z, θ, Z, args...),
        ∇u
    )
end
