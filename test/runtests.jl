using Aqua
using ConstitutiveModels
# using Enzyme
# using ForwardDiff
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
  # P_ad = Tensors.gradient(x -> energy(model, props, x), F)
  # P_enz = pk1_stress(Reverse, model, props, F)
  P_ad = pk1_stress(ADMode, model, props, F)
  # @test P ≈ P_ad
  @test P ≈ P_ad

  # σ = cauchy_stress(model, props, F)
  # # σ_ad = symmetric(J^-1 * dot(Tensors.gradient(x -> energy(model, props, x), F), F'))
  # σ_enz = cauchy_stress(Reverse, model, props, F)

  # # @test σ ≈ σ_ad
  # @test σ ≈ σ_enz
end

function ad_test(
  model::Mod, props::V1,
  motion::Type{Motion}, vals
) where {Mod <: ConstitutiveModel, Motion <: SimpleMotion,
         V1 <: AbstractArray{<:Number, 1}}

  for val in vals
    if motion <: UniaxialStressDisplacementControl
      F = ConstitutiveModels.deformation_gradient(motion, model, props, val, Tensor)
    else
      F = ConstitutiveModels.deformation_gradient(motion, val)
    end
    ad_test(model, props, F)
  end
end

@testset ExtendedTestSet "ConstitutiveModels.jl" begin
  @testset ExtendedTestSet "MechanicalModels" begin
    @testset ExtendedTestSet "HyperelasticModels" begin
      include("hyperelastic_models/test_LinearElastic.jl")
      include("hyperelastic_models/test_NeoHookean.jl")
    end
  end
end

@testset ExtendedTestSet "Aqua" begin
  Aqua.test_all(ConstitutiveModels; ambiguities=false)
end

@testset ExtendedTestSet "JET" begin
  JET.test_package(ConstitutiveModels; target_defined_modules=true)
end