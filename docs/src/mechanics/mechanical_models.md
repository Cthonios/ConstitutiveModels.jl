# Mechanical Models
This set of models are useful for composing multiphysics constitutive models in the presence of mechanical loads as well as by themselves when solving e.g. total-Lagrange or update-Lagrange formulations of solid mechanics.

Balance of Mass

Total Lagrange
``\rho_R = J\rho``

Total Lagrange

``\rho_R\frac{\partial^2\mathbf{u}}{\partial t^2} = \nabla\cdot\mathbf{P} + \rho\mathbf{b}``

Updated Lagrange


## Abstract Interface
```@autodocs
Modules = [ConstitutiveModels]
Order   = [:type, :function]
Pages   = ["MechanicalModels.jl"]
```
```@autodocs
Modules = [ConstitutiveModels]
Order   = [:type, :function]
Pages   = ["HyperelasticModels.jl"]
```