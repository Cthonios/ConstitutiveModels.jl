# Elastic Constants

```@meta
CurrentModule = ConstitutiveModels
```

# Elastic Constants

Many constitutive models for solid mechanics require elastic material properties. While there are numerous ways to characterize isotropic linear elasticity, only two independent elastic constants are required to uniquely define the material response.

`ConstitutiveModels.jl` provides the `ElasticConstants` type to simplify this process. Users may specify any valid pair of isotropic elastic constants, and the remaining constants are computed automatically. This allows constitutive models to work with a consistent internal representation regardless of how the material properties were originally specified.

## Overview

The `ElasticConstants` type stores the five most commonly used elastic constants for isotropic materials:

| Symbol | Name |
|--------|------|
| \(E\) | Young's modulus |
| \(\nu\) | Poisson's ratio |
| \(\kappa\) | Bulk modulus |
| \(\lambda\) | Lamé's first constant |
| \(\mu\) | Shear modulus (Lamé's second constant) |

Internally, every `ElasticConstants` object contains all five values, even though only two independent constants are required during construction.

```julia
ElasticConstants{T}
```

stores

```julia
E
ν
κ
λ
μ
```

where `T` is any subtype of `Number`.

## Constructing Elastic Constants

The primary constructor accepts a dictionary of user-supplied material parameters.

```julia
elastic = ElasticConstants(params)
```

where `params` is a `Dict{String}` containing **exactly two** independent elastic constants.

For example,

```julia
params = Dict(
    "Young's modulus" => 210e9,
    "Poisson's ratio" => 0.30,
)

elastic = ElasticConstants(params)
```

automatically computes

- bulk modulus,
- Lamé's first constant,
- shear modulus,

from the supplied Young's modulus and Poisson's ratio.

## Supported Property Names

The constructor recognizes the following property names.

| Property | Dictionary Key |
|-----------|----------------|
| Young's modulus | `"Young's modulus"` |
| Poisson's ratio | `"Poisson's ratio"` |
| Bulk modulus | `"bulk modulus"` |
| Lamé's first constant | `"Lamé's first constant"` |
| Shear modulus | `"shear modulus"` |

The property names must match these strings exactly.

## Supported Input Combinations

Any valid pair of independent isotropic elastic constants may be supplied.

Currently supported combinations are

| First Constant | Second Constant |
|----------------|-----------------|
| Young's modulus | Poisson's ratio |
| Young's modulus | Bulk modulus |
| Young's modulus | Lamé's first constant |
| Young's modulus | Shear modulus |
| Poisson's ratio | Bulk modulus |
| Poisson's ratio | Lamé's first constant |
| Poisson's ratio | Shear modulus |
| Bulk modulus | Lamé's first constant |
| Bulk modulus | Shear modulus |
| Lamé's first constant | Shear modulus |

These correspond to all ten possible independent pairs of the five isotropic elastic constants.

## Conversion Relations

The constructor automatically applies the appropriate analytical conversion formulas for the supplied pair of elastic constants.

For example, when Young's modulus and Poisson's ratio are provided,

```math
\mu = \frac{E}{2(1+\nu)},
```

```math
\lambda = \frac{E\nu}{(1+\nu)(1-2\nu)},
```

```math
\kappa = \frac{E}{3(1-2\nu)}.
```

Similarly, all other supported input combinations are converted internally using the standard isotropic elasticity identities.

This allows users to specify material data in whichever form is most convenient while ensuring a consistent representation throughout the constitutive model library.

## Accessing Elastic Constants

After construction, each elastic constant is available as a field.

```julia
elastic = ElasticConstants(params)

elastic.E
elastic.ν
elastic.κ
elastic.λ
elastic.μ
```

This avoids repeated conversions inside constitutive models and provides direct access to whichever constant is most convenient for a given formulation.

## Display

Printing an `ElasticConstants` object produces a formatted summary of all computed elastic constants.

For example,

```julia
println(elastic)
```

produces output similar to

```text
Elastic constants:
  bulk modulus          = ...
  Lamé's first constant = ...
  Poisson's ratio       = ...
  shear modulus         = ...
  Young's Modulus       = ...
```

This is useful when verifying material parameters or debugging constitutive models.

## Error Handling

Construction requires **two independent elastic constants**.

If fewer than two constants are supplied, an `ElasticConstantsError` is thrown with a descriptive message indicating the missing information.

For example,

```julia
params = Dict(
    "Young's modulus" => 210e9,
)

ElasticConstants(params)
```

produces an error explaining that a second independent elastic constant is required.

If an unrecognized property name is supplied, the error message also lists the valid property names accepted by the constructor.

This validation helps identify typographical errors and incomplete material definitions during model setup.

## Design Philosophy

Many constitutive models can be formulated using different elastic constants. For example, linear elasticity is often expressed using Lamé parameters, while engineering handbooks typically tabulate Young's modulus and Poisson's ratio. Other applications may naturally use the bulk and shear moduli.

Rather than requiring every constitutive model to perform its own conversions—or forcing users to provide a specific pair of parameters—`ElasticConstants` centralizes these relationships in a single, reusable utility.

This design offers several advantages:

- users may specify whichever pair of elastic constants is most readily available,
- constitutive models always receive a complete and consistent representation,
- analytical conversion formulas are implemented in one location,
- input validation and error reporting are standardized across the package.

By separating elastic parameter conversion from constitutive model implementation, `ConstitutiveModels.jl` provides a flexible interface that accommodates a wide range of engineering workflows while keeping constitutive models simple and maintainable.

```@autodocs
Modules = [ConstitutiveModels]
Order   = [:type, :function]
Pages   = ["ElasticConstants.jl"]
```
