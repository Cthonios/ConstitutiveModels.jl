struct FeFv{
    NP, NS,
    EqStrainEnergyDensity <: AbstractHyperelastic
} <: AbstractMechanicalModel{NP, NS}
    eq_strain_energy_density::EqStrainEnergyDensity
end
