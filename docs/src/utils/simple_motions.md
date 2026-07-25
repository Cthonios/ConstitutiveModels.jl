# Simple motions

```@meta
CurrentModule = ConstitutiveModels
```

# Material Point Simulation and Analytic Motions

One of the primary goals of `ConstitutiveModels.jl` is to make the development and verification of constitutive models as straightforward as possible. Before integrating a constitutive model into a finite element code, it is often desirable to evaluate its response at a **single material point** subjected to prescribed deformation histories.

To support this workflow, `ConstitutiveModels.jl` provides

- a general material point simulator,
- a hierarchy of prescribed kinematic motions,
- constrained loading paths,
- automatic history storage,
- optional temperature and temperature gradient histories.

Together these utilities provide a lightweight "virtual material testing machine" for evaluating constitutive behavior under a variety of loading conditions.

---

# Overview

The simulation framework consists of four primary components.

| Component | Purpose |
|-----------|---------|
| `AbstractMotion` | Represents a prescribed loading path |
| `simulate_material_point` | Advances the constitutive model through time |
| `MaterialHistory` | Stores the complete simulation history |
| Motion implementations | Define common analytical deformation histories |

The general workflow is

```text
Motion
      ↓
Kinematics
      ↓
Constitutive Model
      ↓
Material Response
      ↓
MaterialHistory
```

The same simulation driver can be used for nearly every constitutive model in the package.

---

# Motion Interface

Every loading path derives from

```julia
abstract type AbstractMotion end
```

which represents a time-dependent deformation history.

Two specialized subclasses are provided.

```julia
AbstractSimpleMotion
```

describes motions whose kinematics are prescribed analytically.

Examples include

- uniaxial strain,
- simple shear,
- biaxial stretch,
- pure shear.

These motions require no iterative solution.

---

```julia
AbstractConstrainedMotion
```

represents loading paths in which one or more kinematic quantities must be determined implicitly by satisfying mechanical constraints.

For example,

- uniaxial stress,
- plane stress,
- traction-free lateral boundaries.

These motions generally require a nonlinear solve at every time increment.

---

# Kinematic Quantities

Different constitutive models operate on different kinematic measures.

The simulation framework automatically computes the required quantity.

Depending on the constitutive model, this may be

- displacement gradient

```math
\nabla\mathbf{u},
```

- deformation gradient

```math
\mathbf{F},
```

- velocity gradient

```math
\mathbf{L},
```

- infinitesimal strain

```math
\boldsymbol{\varepsilon}.
```

The helper function

```julia
_kinematics(...)
```

selects the appropriate measure based on the constitutive model.

This allows the same motion definition to be used for

- finite strain hyperelasticity,
- infinitesimal elasticity,
- hypoelasticity,
- elastoplasticity,
- viscoelasticity,

without modification.

---

# Material Point Simulation

The primary user interface is

```julia
simulate_material_point(...)
```

which integrates a constitutive model over a prescribed loading history.

The simulator

1. initializes material properties,
2. initializes internal variables,
3. advances time,
4. evaluates the constitutive response,
5. updates history variables,
6. stores the results.

The result is a vector of `MaterialHistory` objects that completely describe the evolution of the material point.

---

# Hyperelastic Models

For hyperelastic constitutive models, the simulator evaluates constitutive quantities directly from the current deformation.

Each time increment performs

```text
time
      ↓
temperature
      ↓
kinematics
      ↓
constitutive model
      ↓
state update
      ↓
history output
```

Internal state variables are updated after every increment before advancing to the next load step.

---

# Hypoelastic Models

Hypoelastic constitutive models differ slightly because the constitutive response depends on the previous stress state.

Consequently, the simulator additionally stores

```math
\boldsymbol{\sigma}_n
```

and updates

```math
\boldsymbol{\sigma}_{n+1}
```

after each increment.

Apart from this difference, both simulation drivers expose the same interface to the user.

---

# Material History

Simulation results are stored in

```julia
MaterialHistory
```

which records the complete state of the material point at every time step.

Each history entry contains

| Quantity | Description |
|----------|-------------|
| Time | Current simulation time |
| Kinematics | Deformation measure supplied to the constitutive model |
| Material output | Quantity returned by the requested constitutive function |
| State variables | Internal variables after the update |
| Temperature | Current temperature |
| Temperature gradient | Current temperature gradient |

The returned history vector may be postprocessed to generate stress-strain curves, internal variable evolution, energy histories, or other constitutive responses.

---

# Temperature Histories

The simulator optionally accepts prescribed thermal loading.

Temperature is specified by

```julia
temp_func(t)
```

while the temperature gradient is prescribed through

```julia
temp_grad_func(t)
```

For purely mechanical simulations, both default to zero.

Thermomechanical constitutive models can therefore be exercised using exactly the same simulation framework as isothermal models.

---

# Analytical Motion Definitions

`ConstitutiveModels.jl` provides several commonly used analytical loading paths.

These are intended primarily for

- constitutive model verification,
- regression testing,
- parameter calibration,
- visualization,
- comparison with analytical solutions.

---

# Uniaxial Strain

The uniaxial strain motion prescribes

```math
\mathbf{F}
=
\begin{bmatrix}
\lambda(t) & 0 & 0\\
0 & 1 & 0\\
0 & 0 & 1
\end{bmatrix}.
```

Only the stretch history

```math
\lambda(t)
```

must be provided.

This loading path is commonly used for

- tensile testing,
- compression,
- elastoplastic verification.

Example

```julia
motion = UniaxialStrain(t -> 1 + 0.2t)
```

---

# Biaxial Strain

Biaxial strain independently stretches two principal directions.

```math
\mathbf{F}
=
\begin{bmatrix}
\lambda_1(t) & 0 & 0\\
0 & \lambda_2(t) & 0\\
0 & 0 & 1
\end{bmatrix}.
```

This loading path is useful for

- rubber elasticity,
- membrane materials,
- anisotropic constitutive models.

Example

```julia
motion = BiaxialStrain(
    t -> 1 + 0.2t,
    t -> 1 + 0.1t
)
```

---

# Isochoric Uniaxial Stress

This motion preserves volume while stretching one principal direction.

The deformation gradient is

```math
\mathbf{F}
=
\begin{bmatrix}
\lambda & 0 & 0\\
0 & \lambda^{-1/2} & 0\\
0 & 0 & \lambda^{-1/2}
\end{bmatrix},
```

which satisfies

```math
\det(\mathbf{F})=1.
```

Isochoric loading is frequently used when testing nearly incompressible constitutive models.

---

# Pure Shear

Pure shear is represented by

```math
\mathbf{F}
=
\frac12
\begin{bmatrix}
\lambda+\lambda^{-1} &
\lambda-\lambda^{-1} &
0\\
\lambda-\lambda^{-1} &
\lambda+\lambda^{-1} &
0\\
0&0&2
\end{bmatrix}.
```

Unlike simple shear, this loading path consists entirely of principal stretches with no rigid-body rotation.

It is widely used in nonlinear elasticity.

---

# Simple Shear

Simple shear prescribes

```math
\mathbf{F}
=
\begin{bmatrix}
1&\gamma(t)&0\\
0&1&0\\
0&0&1
\end{bmatrix}.
```

where

```math
\gamma(t)
```

is the prescribed engineering shear.

Simple shear is one of the standard benchmark problems for

- finite deformation elasticity,
- viscoelasticity,
- crystal plasticity,
- large-strain plasticity.

Example

```julia
motion = SimpleShear(t -> 0.5t)
```

---

# Velocity Gradient

Every analytical motion also provides its associated velocity gradient

```math
\mathbf{L}
=
\dot{\mathbf{F}}\mathbf{F}^{-1},
```

computed automatically through automatic differentiation using `ForwardDiff.jl`.

Consequently, the user only specifies the deformation history, while the corresponding velocity gradient is generated automatically.

This allows the same motion to drive both

- finite deformation formulations using the deformation gradient, and
- rate-based constitutive models requiring the velocity gradient.

---

# Constrained Motions

Not every experimental loading condition can be prescribed analytically.

For example, a uniaxial tensile test is typically performed under

```math
\sigma_{22}=\sigma_{33}=0,
```

rather than prescribing the lateral stretches.

The class

```julia
AbstractConstrainedMotion
```

supports these situations.

---

# Uniaxial Stress Displacement Control

The current implementation provides

```julia
UniaxialStressDisplacementControl
```

which prescribes the axial displacement while solving for the unknown lateral strains.

At each time increment, a nonlinear system is solved such that

```math
\sigma_{22}=0,
```

```math
\sigma_{33}=0.
```

The unknown lateral deformation is determined using a Newton iteration based on automatic differentiation.

This makes the loading path applicable to arbitrary nonlinear constitutive models without requiring analytical tangent expressions.

---

# Newton Solver

Constrained motions rely on an internal Newton solver.

Given a nonlinear residual

```math
\mathbf{r}(\mathbf{x}),
```

the update is computed as

```math
\Delta\mathbf{x}
=
-\mathbf{J}^{-1}\mathbf{r},
```

where

```math
\mathbf{J}
=
\frac{\partial\mathbf{r}}{\partial\mathbf{x}}
```

is evaluated automatically using `ForwardDiff.jl`.

Iterations continue until both the residual and solution updates satisfy the prescribed convergence tolerances.

Although this solver is primarily an implementation detail, it enables constrained loading paths to remain independent of any particular constitutive model.

---

# Design Philosophy

The material point simulator is intended to provide a common testing environment for every constitutive model in `ConstitutiveModels.jl`.

By separating

- constitutive behavior,
- kinematic loading paths,
- state evolution,
- temperature histories,
- and result storage,

the same simulation framework can be reused across a wide range of material models with minimal code duplication.

This modular design also makes it straightforward to add new analytical motions, constrained loading conditions, or constitutive models without modifying the simulation driver itself, providing a flexible foundation for the multiphysics capabilities of `ConstitutiveModels.jl`.

```@autodocs
Modules = [ConstitutiveModels]
Order   = [:type, :function]
Pages   = ["SimpleMotions.jl"]
```
