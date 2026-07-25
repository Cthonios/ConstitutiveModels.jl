# Material Symmetry

```@meta
CurrentModule = ConstitutiveModels
```

# Material Symmetry

The `MaterialSymmetry` interface in `ConstitutiveModels.jl` provides a unified mechanism for representing constitutive parameters with different tensor symmetries. Rather than requiring every constitutive model to explicitly construct second- or fourth-order tensors, the symmetry object is responsible for converting a compact list of material properties into the appropriate tensor representation.

This separation provides several advantages:

- Constitutive models remain independent of the chosen material symmetry.
- Material parameters can be stored compactly.
- New symmetry classes can be introduced without modifying existing constitutive models.
- The same constitutive model can be used with isotropic, cubic, orthotropic, transversely isotropic, or fully anisotropic material parameters.

The current implementation provides support for isotropic second- and fourth-order tensors.

## Interface

All material symmetry objects inherit from

```julia
abstract type AbstractMaterialSymmetry{O} end
```

where `O` denotes the tensor order represented by the symmetry object.

Every material symmetry implementation is expected to provide three methods.

```julia
as_tensor(::AbstractMaterialSymmetry, props)
```

Constructs the tensor representation from a compact vector of material properties.

```julia
initialize_props(::AbstractMaterialSymmetry, inputs, keys)
```

Initializes the material property vector from user-provided inputs.

```julia
num_properties(::AbstractMaterialSymmetry)
```

Returns the number of independent material properties required by the symmetry.

These functions define the complete interface required by constitutive models.

## Tensor Order

The type parameter `O` specifies the tensor order.

| Order | Meaning | Typical Use |
|-------:|---------|-------------|
| `2` | Second-order symmetric tensor | Thermal conductivity, dielectric permittivity, magnetic permeability, diffusion tensors |
| `4` | Fourth-order symmetric tensor | Elastic stiffness, compliance, tangent moduli |

The tensor order is encoded in the type rather than stored as runtime information, allowing Julia to specialize and optimize tensor construction.

## Isotropic Symmetry

An isotropic material possesses identical properties in every direction. Consequently, only a small number of independent material parameters are required.

### Second-Order Isotropy

A second-order isotropic tensor has the form

```math
\mathbf{A} = a\mathbf{I},
```

where

- \(a\) is the single material parameter,
- \(\mathbf{I}\) is the second-order identity tensor.

Accordingly, only one material property is required.

```julia
Isotropy{2}()
```

requires

```julia
num_properties(Isotropy{2}()) == 1
```

Internally,

```julia
as_tensor(Isotropy{2}(), props)
```

constructs

```math
A_{ij}=a\,\delta_{ij}.
```

This representation is commonly used for

- isotropic thermal conductivity,
- isotropic electrical conductivity,
- isotropic magnetic permeability,
- isotropic diffusivity.

---

### Fourth-Order Isotropy

A fourth-order isotropic tensor is represented as

```math
C_{ijkl}
=
\lambda\,\delta_{ij}\delta_{kl}
+
\mu
\left(
\delta_{ik}\delta_{jl}
+
\delta_{il}\delta_{jk}
\right),
```

where

- \(\lambda\) is the first Lamé parameter,
- \(\mu\) is the shear modulus,
- \(\delta_{ij}\) is the Kronecker delta.

Therefore, two independent material constants are required.

```julia
Isotropy{4}()
```

requires

```julia
num_properties(Isotropy{4}()) == 2
```

and

```julia
as_tensor(Isotropy{4}(), props)
```

constructs the fourth-order symmetric tensor from the supplied Lamé parameters.

This representation is used extensively for

- linear elasticity,
- tangent stiffness operators,
- isotropic constitutive tangents.

## Property Initialization

Material properties are initialized through

```julia
initialize_props(symmetry, inputs, keys)
```

where

- `inputs` contains the user-supplied material data,
- `keys` specifies which properties should be extracted.

For example,

```julia
symmetry = Isotropy{4}()

props = initialize_props(
    symmetry,
    inputs,
    [:λ, :μ],
)
```

produces a compact property vector

```julia
[λ, μ]
```

which can later be converted into a tensor via

```julia
C = as_tensor(symmetry, props)
```

Separating initialization from tensor construction allows constitutive models to remain independent of the storage format of material parameters.

## Extending the Interface

New material symmetry classes can be implemented by inheriting from `AbstractMaterialSymmetry{O}` and defining the three required interface functions.

For example,

```julia
struct Cubic{O} <: AbstractMaterialSymmetry{O}
end
```

would provide specialized implementations of

```julia
initialize_props(::Cubic, ...)
```

```julia
as_tensor(::Cubic, ...)
```

```julia
num_properties(::Cubic)
```

without requiring any changes to the constitutive models themselves.

This design makes it straightforward to add support for additional symmetry classes such as

- cubic,
- transversely isotropic,
- orthotropic,
- monoclinic,
- triclinic,

or application-specific material symmetries.

## Design Philosophy

The material symmetry abstraction is intended to separate **material symmetry** from **constitutive behavior**.

A constitutive model describes the physical relationship between stress, strain, temperature, electric field, magnetic field, or other state variables, while the material symmetry determines how the associated constitutive tensors are represented.

This modular design enables the same constitutive model implementation to be reused with multiple symmetry classes, reducing code duplication and making it easier to extend `ConstitutiveModels.jl` to new classes of multiphysics material models.
