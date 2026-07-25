# Interface

```@meta
CurrentModule = ConstitutiveModels
```

# Constitutive Model Interface

`ConstitutiveModels.jl` is built around a small set of composable interfaces that separate **material behavior**, **kinematics**, **material symmetry**, and **constitutive modules** into independent components.

Rather than requiring every constitutive model to implement a large collection of methods, the package provides sensible defaults wherever possible. Many constitutive quantities—including stresses, tangent moduli, entropy, and heat capacity—can be generated automatically using automatic differentiation.

This design minimizes boilerplate while allowing constitutive models to remain concise, extensible, and physically transparent.

---

# Package Architecture

At the highest level, every object participating in the constitutive framework derives from

```julia
AbstractConstitutive
```

The hierarchy is

```text
AbstractConstitutive
│
├── AbstractConstitutiveModel
│   ├── AbstractHyperelasticModel
│   ├── AbstractHypoelasticModel
│   └── AbstractLinearElasticModel
│
├── AbstractConstitutiveModule
│
└── AbstractMaterialSymmetry
```

Each level introduces progressively more specialized functionality while remaining compatible with the rest of the framework.

---

# AbstractConstitutive

Every constitutive object derives from

```julia
abstract type AbstractConstitutive end
```

This common ancestor defines the minimum interface expected throughout the package.

Every constitutive object should implement

```julia
initialize_props(...)
```

```julia
num_properties(...)
```

```julia
num_state_variables(...)
```

and may optionally implement

```julia
state_variable_names(...)
```

The package provides sensible default implementations wherever possible.

---

# Material Properties

Every constitutive object is responsible for converting user-supplied material parameters into a compact property vector.

This is accomplished through

```julia
initialize_props(model, inputs)
```

which transforms a dictionary of user inputs into an ordered vector suitable for efficient numerical evaluation.

The number of required material parameters is given by

```julia
num_properties(model)
```

This interface allows the package to validate input data before constitutive calculations begin.

---

# Internal State Variables

Many constitutive models contain internal variables describing irreversible processes such as

- plastic strain,
- hardening,
- damage,
- viscoelastic internal variables,
- phase fractions.

The required number of state variables is reported by

```julia
num_state_variables(model)
```

By default,

```julia
initialize_state(model)
```

returns

```julia
zeros(num_state_variables(model))
```

which is appropriate for many constitutive models.

Models requiring more sophisticated initialization may override this method.

---

# State Variable Names

To improve postprocessing and visualization, constitutive models may assign human-readable names to their state variables.

```julia
state_variable_names(model)
```

returns a vector of strings describing each stored variable.

The default implementation generates

```text
state_1
state_2
state_3
...
```

although constitutive models are encouraged to provide more descriptive names such as

```text
Equivalent Plastic Strain
Backstress
Damage
```

This information is particularly useful when exporting material histories from simulations.

---

# Constitutive Models

Actual constitutive laws derive from

```julia
AbstractConstitutiveModel
```

These objects describe the physical relationship between deformation, temperature, internal variables, and stress.

All constitutive models share a common calling convention

```julia
method(
    model,
    props,
    Z_old,
    Z_new,
    Δt,
    ...
)
```

where

| Argument | Description |
|----------|-------------|
| `model` | Constitutive model |
| `props` | Material properties |
| `Z_old` | Previous state variables |
| `Z_new` | Updated state variables |
| `Δt` | Time increment |

Additional arguments depend on the model type.

---

# Constitutive Model Categories

The package distinguishes constitutive models according to the kinematic quantity they require.

## Hyperelastic Models

```julia
AbstractHyperelasticModel
```

receive the displacement gradient

```math
\nabla\mathbf{u}.
```

These models are typically defined through a Helmholtz free energy function

```math
\psi(\nabla\mathbf{u},\theta,\mathbf{Z}).
```

Examples include

- nonlinear elasticity,
- hyperelasticity,
- finite-strain thermoelasticity.

---

## Hypoelastic Models

```julia
AbstractHypoelasticModel
```

receive

- the velocity gradient,
- the previous Cauchy stress.

These models evolve stress incrementally and are generally used for objective stress-rate formulations.

---

## Linear Elastic Models

```julia
AbstractLinearElasticModel
```

operate directly on the infinitesimal strain tensor

```math
\boldsymbol{\varepsilon}.
```

They share much of the hyperelastic interface but assume small deformations.

---

# Automatic Differentiation

One of the central design principles of `ConstitutiveModels.jl` is that many constitutive quantities can be computed automatically from the Helmholtz free energy.

If a constitutive model implements

```julia
helmholtz_free_energy(...)
```

the package automatically provides default implementations of

- First Piola stress
- Cauchy stress
- Material tangent
- Entropy
- Heat capacity
- Thermomechanical coupling tensors

using automatic differentiation through `Tensors.jl`.

For example,

```math
\mathbf{P}
=
\frac{\partial\psi}{\partial\mathbf{F}},
```

is computed automatically by

```julia
pk1_stress(...)
```

Similarly,

```math
\eta
=
-
\frac{\partial\psi}{\partial\theta},
```

defines the entropy, while

```math
c
=
\theta
\frac{\partial\eta}{\partial\theta}
```

gives the heat capacity.

This approach greatly reduces the amount of constitutive code that must be written and minimizes opportunities for implementation errors.

---

# Density

Every constitutive model is expected to provide a reference density.

By default,

```julia
density(...)
```

returns the first entry of the material property vector.

Consequently, the first material property is conventionally the Lagrangian density.

Models with more sophisticated density evolution may override this method.

---

# Constitutive Modules

Many constitutive models are naturally composed from smaller physical building blocks.

These reusable components derive from

```julia
AbstractConstitutiveModule
```

Examples include

- hardening laws,
- yield surfaces,
- thermal expansion,
- damage models,
- equation-of-state models.

Modules possess their own

- material properties,
- state variables,
- initialization routines,

allowing them to be reused across multiple constitutive models.

The helper

```julia
module_props(...)
```

extracts the appropriate slice of the global material property vector corresponding to a particular module.

This enables modular constitutive models without sacrificing storage efficiency.

---

# Material Symmetry

Material symmetry objects are independent of constitutive behavior.

Rather than storing full constitutive tensors, a symmetry object converts compact parameter vectors into tensor representations.

Examples include

- isotropic,
- cubic,
- orthotropic,
- transversely isotropic,

and future symmetry classes.

Separating symmetry from constitutive behavior allows the same constitutive model to be used with multiple material classes.

---

# Kinematics

Different constitutive models require different measures of deformation.

The package introduces lightweight marker types

- `DisplacementGradient`
- `VelocityGradient`
- `LinearStrain`

to identify the expected kinematics.

The function

```julia
kinematics(model)
```

returns the appropriate measure automatically based on the model type.

Simulation drivers and finite element codes therefore do not need to know the specific requirements of individual constitutive models.

---

# Temperature Gradient Dependencies

Most constitutive quantities depend only on

- deformation,
- temperature,
- internal variables.

A small number additionally require the temperature gradient, such as

- heat flux,
- dissipation.

The helper

```julia
requires_temperature_gradient(...)
```

allows algorithms to determine automatically whether the additional argument is needed, avoiding unnecessary computation for purely mechanical constitutive quantities.

---

# Design Philosophy

The guiding philosophy of `ConstitutiveModels.jl` is to separate independent physical concepts into reusable interfaces.

- **Constitutive models** describe material behavior.
- **Constitutive modules** describe reusable physical mechanisms.
- **Material symmetries** describe tensor representations.
- **Kinematic types** describe the required deformation measures.
- **Automatic differentiation** derives many constitutive quantities directly from thermodynamic potentials.

This modular architecture allows complex multiphysics constitutive models to be assembled from small, composable components while minimizing duplicated code and ensuring consistency throughout the package.

```@autodocs
Modules = [ConstitutiveModels]
Order   = [:type, :function]
Pages   = ["Interface.jl"]
```