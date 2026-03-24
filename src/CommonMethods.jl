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
    props, Δt, Z_old, Z_new,
    ∇u, θ, args...
)
    F = ∇u + one(∇u)
    J = det(F)
    P = pk1_stress(model, props, Δt, Z_old, Z_new, ∇u, θ, args...)
    σ = (1 / J) * symmetric(P ⋅ transpose(F))
    return σ
end

function entropy(
    model::ModelsWithThermal,
    props, Δt, Z_old, Z_new,
    ∇u, θ, args...
)
    return entropy(
        ForwardDiffAD(),
        model, props, Δt,
        ∇u, θ, Z_old, Z_new
    )
end

function heat_capacity(
    model::ModelsWithThermal,
    props, Δt, Z_old, Z_new,
    ∇u, θ, args...
)
    return heat_capacity(
        ForwardDiffAD(),
        model, props, Δt, Z_old, Z_new,
        ∇u, θ, args...
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
    props, Δt, Z_old, Z_new,
    ∇u, θ, args...
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
    props, Δt, Z_old, Z_new,
    ∇u, θ, args...
)
    return material_tangent(
        ForwardDiffAD(),
        model, props, Δt, Z_old, Z_new,
        ∇u, θ, args...
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
    props, Δt, Z_old, Z_new,
    ∇u, θ, args...
)
    return pk1_stress(
        ForwardDiffAD(),
        model, props, Δt, Z_old, Z_new,
        ∇u, θ, args...
    )
end

# ForwardDiff wrappers

function entropy(
    ::ForwardDiffAD,
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    η = -Tensors.gradient(z -> 
        helmholtz_free_energy(model, props, Δt, ∇u, z, Z_old, Z_new),
        θ
    )
    return η
end

function heat_capacity(
    ::ForwardDiffAD,
    model::ModelsWithThermal,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    c = -θ * Tensors.hessian(z -> 
        helmholtz_free_energy(model, props, Δt, ∇u, z, Z_old, Z_new),
        θ
    )
    return c
end

function material_tangent(
    ::ForwardDiffAD,
    model::ModelsWithMechanics,
    props, Δt,
    ∇u, θ, Z_old, Z_new
)
    return Tensors.gradient(z -> 
        pk1_stress(ForwardDiffAD(), model, props, Δt, z, θ, Z_old, Z_new),
        ∇u
    )
end

function pk1_stress(
    ::ForwardDiffAD,
    model::ModelsWithMechanics,
    props, Δt, Z_old, Z_new,
    ∇u, θ, args...
)
    return Tensors.gradient(z -> 
        helmholtz_free_energy(model, props, Δt, Z_old, Z_new, z, θ, args...),
        ∇u
    )
end
