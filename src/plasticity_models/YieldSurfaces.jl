abstract type YieldSurface{NProps, NStateVars} <: PlasticityModel{NProps, NStateVars} end

# 1 prop for yield stress
# 1 props for radius
struct J2YieldSurface <: YieldSurface{1, 1}
end

function J2YieldSurface(inputs::D) where D <: Dict
  @assert "yield stress" in keys(inputs)
  yield_stress = inputs["yield stress"]
  return J2YieldSurface(), @SVector [yield_stress]
end