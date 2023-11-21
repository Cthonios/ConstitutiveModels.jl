using Aqua
using ConstitutiveModels
using Enzyme
using ForwardDiff
using JET
using Tensors
using Test
using TestSetExtensions

# function props_test(inputs::Dict{String, <:Number}, props::V) where V <: AbstractArray{<:Number, 1}

# end

function ad_test(model::Mod, props::V, F::Tensor{2, 3, <:Number, 9}) where {Mod <: ConstitutiveModel, V <: AbstractArray{<:Number, 1}}
  C = tdot(F)
  J = det(F)
  
  P = pk1_stress(model, props, F)
  P_ad = Tensors.gradient(x -> energy(model, props, x), F)
  P_enz = pk1_stress(Reverse, model, props, F)
  @test P ≈ P_ad
  @test P ≈ P_enz

  S = pk2_stress(model, props, C)
  S_ad = 2. * Tensors.gradient(x -> energy(model, props, x), C)
  S_enz = pk2_stress(Reverse, model, props, C)
  @test S ≈ S_ad
  @test S ≈ S_enz

  σ = cauchy_stress(model, props, F)
  σ_ad = J^-1 * dot(Tensors.gradient(x -> energy(model, props, x), F), F')
  σ_enz = cauchy_stress(Reverse, model, props, F)
  @test σ ≈ σ_ad
  @test σ ≈ σ_enz
end

function ad_test(
  model::Mod, props::V1,
  motion::Type{Motion}, vals
) where {Mod <: ConstitutiveModel, Motion <: SimpleMotion,
         V1 <: AbstractArray{<:Number, 1}}

  for val in vals
    F = deformation_gradient(motion, val)
    ad_test(model, props, F)
  end
end

@testset ExtendedTestSet "ConstitutiveModels.jl" begin
  @testset ExtendedTestSet "MechanicalModels" begin
    include("mechanical_models/test_NeoHookean.jl")
  end
end

@testset ExtendedTestSet "Aqua" begin
  # Aqua.test_all(ConstitutiveModels)
  Aqua.test_ambiguities(ConstitutiveModels)
  Aqua.test_unbound_args(ConstitutiveModels)
  Aqua.test_undefined_exports(ConstitutiveModels)
  Aqua.test_piracies(ConstitutiveModels)
  Aqua.test_project_extras(ConstitutiveModels)
  Aqua.test_stale_deps(ConstitutiveModels)
  Aqua.test_deps_compat(ConstitutiveModels)
  # Aqua.test_project_toml_formatting(ConstitutiveModels)

end

# @testset ExtendedTestSet "JET" begin
#   JET.test_package(ConstitutiveModels; target_defined_module=true)
# end