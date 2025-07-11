"""
$(TYPEDEF)
"""
abstract type AbstractThermalModel{NP, NS} <: AbstractConstitutiveModel{NP, NS}
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
``\\mathbf{q}``  - heat flux vector\n
``\\mathcal{Z}`` - new state variable array (different from input ```Z```)
$(TYPEDSIGNATURES)
"""
function heat_flux(
    ::AbstractThermalModel{NP, NS},
    props,
    Δt,
    ∇θ,
    θ,
    Z,
    args...
) where {NP, NS}
    @assert false "needs to be implemented for new model"
end

include("FouriersLaw.jl")
