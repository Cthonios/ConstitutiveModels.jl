abstract type HardeningModel{NProps, NStateVars} <: PlasticityModel{NProps, NStateVars} end
abstract type IsotropicHardeningModel{NProps, NStateVars} <: HardeningModel{NProps, NStateVars} end


struct NoIsotropicHardening <: IsotropicHardeningModel{0, 0}
end

function NoIsotropicHardening(::D) where D <: Dict
  return NoIsotropicHardening(), SVector{0, Float64}(), SVector{0, Float64}()
end

# hardening(::NoIsotropicHardening, )
radius(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray = sqrt(2. / 3.) * props[1]
slope(::NoIsotropicHardening, props::V, eqps) where V <: AbstractArray  = 0.0