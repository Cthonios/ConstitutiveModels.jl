struct MaxwellElement{
    NP, NS,
    H <: AbstractHyperelastic,
    V <: AbstractViscosity
} <: AbstractViscoelasticity{NP, NS, H, V}
    spring::H
    dashpot::V
end

