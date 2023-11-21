struct LinearElastoPlasticity{NProps, Yield <: YieldSurface, Hard <: HardeningModel} <: PlasticityModel{NProps, 7}
  elastic_model::LinearElastic
  yield_surface::Yield
  hardening_model::Hard
end

function LinearElastoPlasticity(inputs::D) where D <: Dict
  elastic_model, elastic_props     = LinearElastic(inputs)
  yield_surface, yield_props       = J2YieldSurface(inputs)
  hardening_model, hardening_props = NoIsotropicHardening(inputs)
  props = vcat(elastic_props, yield_props, hardening_props)
  return LinearElastoPlasticity{3, typeof(yield_surface), typeof(hardening_model)}(
    elastic_model, yield_surface, hardening_model
  ), props
end

function cauchy_stress(model::LinearElastoPlasticity, props, F)
  @show "here"
  σ_el = cauchy_stress(model.elastic_model, props, F)
end

