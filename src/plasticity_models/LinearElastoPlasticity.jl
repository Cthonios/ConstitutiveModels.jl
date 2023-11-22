struct LinearElastoPlasticity{NProps, Yield <: YieldSurface} <: PlasticityModel{NProps, 7}
  elastic_model::LinearElastic
  yield_surface::Yield
end

function LinearElastoPlasticity(inputs::D) where D <: Dict
  elastic_model, elastic_props, elastic_state       = LinearElastic(inputs)
  yield_surface, yield_props, yield_state           = J2YieldSurface(inputs)
  props = vcat(elastic_props, yield_props)
  state = vcat(elastic_state, yield_state)
  return LinearElastoPlasticity{number_of_properties(elastic_model) + number_of_properties(yield_surface), 
                                typeof(yield_surface)}(
    elastic_model, yield_surface
  ), props, state
end

function yield_surface_properties(model::LinearElastoPlasticity, props)
  n_props = number_of_properties(model) - number_of_properties(model.elastic_model)
  return @views SVector{n_props, eltype(props)}(props[3:end])
end

function cauchy_stress(model::LinearElastoPlasticity, props, F, state)
  # extract props
  μ           = props[2]
  yield_props = yield_surface_properties(model, props)

  # kinematics
  I = one(SymmetricTensor{2, 3, eltype(F), 6})
  ∇u = F - I
  ε = symmetric(∇u)

  # unpack state variables
  ε_p_old = @views fromvoigt(SymmetricTensor{2, 3, eltype(ε), 6}, state[1:6])
  α_old   = state[7]

  # calculate elastic trial stress
  ε_e_tr = ε - ε_p_old
  σ_e_tr = cauchy_stress(model.elastic_model, elastic_properties(props), ε_e_tr, state)

  # calculate increment in plastic flow
  Δγ = update(model.yield_surface, yield_props, μ, σ_e_tr, α_old)

  # update stuff
  if Δγ > 0.0
    N = dev(σ_e_tr) / norm(dev(σ_e_tr))
    σ = σ_e_tr - 2. * μ * Δγ * N
    ε_p_new = ε_p_old + Δγ * N
    α_new   = α_old   + sqrt(2. / 3.) * Δγ
  else
    σ = σ_e_tr
    ε_p_new = ε_p_old
    α_new   = α_old
  end

  return σ, vcat(tovoigt(SVector, ε_p_new), α_new)
end

function pk1_stress(model::LinearElastoPlasticity, props, F, state)
  J = det(F)
  σ, state = cauchy_stress(model, props, F, state)
  P = J * dot(σ, inv(F)')
  return P
end

