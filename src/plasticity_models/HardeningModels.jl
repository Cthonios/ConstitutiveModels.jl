abstract type HardeningModel{NProps, NStateVars} <: PlasticityModel{NProps, NStateVars} end
abstract type IsotropicHardeningModel{NProps, NStateVars} <: HardeningModel{NProps, NStateVars} end


struct NoIsotropicHardening <: IsotropicHardeningModel{1, 1}
end

function NoIsotropicHardening(::D) where D <: Dict
  return NoIsotropicHardening(), @SVector []
end

initialize_state(::NoIsotropicHardening) = SVector{1, Float64}(0.0)

# function hardening(::NoIsotropicHardening, )