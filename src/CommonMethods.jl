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
    ∇u, θ, Z_old, Z_new,
    args...
)
    F = ∇u + one(∇u)
    J = det(F)
    P = pk1_stress(model, props, Δt, ∇u, θ, Z_old, Z_new, args...)
    # σ = (1. / J) * dot(P, F')
    σ = (1 / J) * P ⋅ transpose(F)
    return σ
end

function entropy(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    return entropy(
        model, props, Δt,
        ∇u, θ, Z_old, Z_new,
        ForwardDiffAD()
    )
end

function heat_capacity(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    return heat_capacity(
        model, props, Δt,
        ∇u, θ, Z_old, Z_new,
        ForwardDiffAD()
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
    ∇u, θ, Z_old, Z_new
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
    ∇u, θ, Z_old, Z_new
)
    return material_tangent(
        model, props, Δt,
        ∇u, θ, Z_old, Z_new,
        ForwardDiffAD()
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
    ∇u, θ, Z_old, Z_new,
)
    return pk1_stress(
        model, props, Δt,
        ∇u, θ, Z_old, Z_new,
        ForwardDiffAD()
    )
end

# ForwardDiff wrappers

function entropy(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    η = -Tensors.gradient(z -> 
        helmholtz_free_energy(model, props, Δt, ∇u, z, Z_old, Z_new),
        θ
    )
    return η
end

function heat_capacity(
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    c = -θ * Tensors.hessian(z -> 
        helmholtz_free_energy(model, props, Δt, ∇u, z, Z_old, Z_new),
        θ
    )
    return c
end

function material_tangent(
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    return Tensors.gradient(z -> 
        pk1_stress(model, props, Δt, z, θ, Z_old, Z_new, ForwardDiffAD()),
        ∇u
    )
end

function pk1_stress(
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z_old, Z_new,
    ::ForwardDiffAD
)
    return Tensors.gradient(z -> 
        helmholtz_free_energy(model, props, Δt, z, θ, Z_old, Z_new),
        ∇u
    )
end
