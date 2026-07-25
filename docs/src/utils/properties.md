# Material Property Utilities

```@meta
CurrentModule = ConstitutiveModels
```

# Material Property Utilities

`ConstitutiveModels.jl` stores material parameters as dictionaries of user-supplied inputs. The property utility functions provide a consistent interface for retrieving these parameters while performing validation, generating informative error messages, and supporting stochastic material properties.

Rather than allowing constitutive models to directly access the input dictionary, all material parameters should be retrieved through `get_property`. This centralizes input validation and ensures consistent behavior throughout the package.

The current implementation supports

- deterministic scalar material properties,
- randomly sampled properties from prescribed probability distributions,
- informative exceptions for missing or malformed inputs.

## Overview

The primary interface is

```julia
get_property(inputs, key)
```

which retrieves the requested property from a dictionary of material parameters.

Depending on the stored value, the function may

- return the scalar value directly,
- generate a random sample from a supported probability distribution,
- throw an informative exception if the input is invalid.

This abstraction allows constitutive models to remain agnostic to how material properties are represented.

## Scalar Properties

The most common use case is a deterministic scalar material parameter.

For example,

```julia
inputs = Dict(
    "Young's Modulus" => 210e9,
    "Poisson Ratio"   => 0.30,
)
```

Retrieving a property is straightforward.

```julia
E = get_property(inputs, "Young's Modulus")
ν = get_property(inputs, "Poisson Ratio")
```

Since both entries are numeric, the stored values are returned directly.

This is the recommended format for deterministic constitutive parameters.

## Random Material Properties

Many applications require uncertainty in constitutive parameters, such as stochastic finite element methods, uncertainty quantification, or Monte Carlo simulations.

To support these workflows, a property may instead be specified as a nested dictionary describing a probability distribution.

Currently, the supported format is a normal (Gaussian) distribution.

```julia
inputs = Dict(
    "Yield Stress" => Dict(
        "mean"     => 450.0,
        "variance" => 25.0,
    )
)
```

Calling

```julia
σy = get_property(inputs, "Yield Stress")
```

constructs a normal distribution

```math
X \sim \mathcal{N}(\mu,\sigma^2)
```

using the supplied mean and variance and returns a randomly sampled realization.

Each call to `get_property` produces a new sample, making this interface convenient for stochastic simulations.

## Missing Properties

If a requested property is not present, `get_property` throws a `PropertyNotFoundError`.

For example,

```julia
inputs = Dict(
    "Young's Modulus" => 210e9,
)

ν = get_property(inputs, "Poisson Ratio")
```

results in an error similar to

```text
Could not find property in inputs with key Poisson Ratio.
Properties passed by user are:
  Young's Modulus
```

Listing the available property names makes typographical mistakes and missing inputs easier to identify.

## Invalid Nested Inputs

Nested dictionaries are interpreted as structured property descriptions.

If the nested dictionary does not match one of the supported formats, a `BadNestedPropertyError` is thrown.

For example,

```julia
inputs = Dict(
    "Yield Stress" => Dict(
        "average" => 450.0,
    )
)
```

will produce an error explaining that the nested input format is not currently supported and will display the expected structure.

At present, the supported nested format is

```julia
Dict(
    "mean"     => <Real>,
    "variance" => <Real>,
)
```

Future versions of `ConstitutiveModels.jl` may extend this interface to support additional probability distributions or more sophisticated parameter descriptions.

## Default Values

`get_property` accepts an optional default value.

```julia
ρ = get_property(inputs, "Density", 7850.0)
```

If the property exists, its value is returned.

If it is missing, the current implementation raises a `PropertyNotFoundError` before the default is used. This behavior intentionally enforces that required constitutive parameters must be explicitly provided by the user.

## Error Types

Two custom exception types are defined by the property utilities.

### `PropertyNotFoundError`

Thrown when a requested property cannot be found in the input dictionary.

The exception reports

- the missing property name,
- all available property names supplied by the user.

This provides considerably more useful diagnostics than Julia's standard dictionary lookup errors.

### `BadNestedPropertyError`

Thrown when a nested dictionary does not match one of the recognized input formats.

The exception reports

- the offending property,
- the malformed nested dictionary,
- the currently supported nested input syntax.

Providing dedicated exception types makes it easier to diagnose user input errors during model setup.

## Design Philosophy

The property utility layer separates **constitutive models** from **input parsing**.

Constitutive models should only request the material parameters they require through `get_property`, without needing to know whether those parameters originate from deterministic values, stochastic distributions, or future input formats.

Centralizing property retrieval provides several benefits:

- consistent validation across the package,
- informative and user-friendly error messages,
- support for stochastic material properties,
- a single location for extending input formats in future releases.

As `ConstitutiveModels.jl` evolves, additional property representations—such as log-normal distributions, uniform distributions, correlated random variables, tabulated data, or parameterized functions—can be introduced by extending `get_property` while leaving constitutive models unchanged.