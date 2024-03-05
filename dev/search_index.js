var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = ConstitutiveModels","category":"page"},{"location":"#ConstitutiveModels","page":"Home","title":"ConstitutiveModels","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ConstitutiveModels.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [ConstitutiveModels]","category":"page"},{"location":"#ConstitutiveModels.AbstractMotion","page":"Home","title":"ConstitutiveModels.AbstractMotion","text":"\n\n\n\n","category":"type"},{"location":"#ConstitutiveModels.IsochoricUniaxialStress","page":"Home","title":"ConstitutiveModels.IsochoricUniaxialStress","text":"Provides an analytic motion for uniaxial stress assuming perfect incompressibility.\n\nThis is\n\nmathbfF = beginbmatrix\nlambda  0                         0 \n0        frac1sqrtlambda  0 \n0        0                         frac1sqrtlambda\nendbmatrix\n\n\n\n\n\n","category":"type"},{"location":"#ConstitutiveModels.SimpleMotion","page":"Home","title":"ConstitutiveModels.SimpleMotion","text":"\n\n\n\n","category":"type"},{"location":"#ConstitutiveModels.SimpleShear","page":"Home","title":"ConstitutiveModels.SimpleShear","text":"Provides an analytic motion for simple shear.\n\nThis is\n\nmathbfF = beginbmatrix\n1  gamma  0 \n0  1       0 \n0  0       1\nendbmatrix\n\n\n\n\n\n","category":"type"},{"location":"#ConstitutiveModels.UniaxialStrain","page":"Home","title":"ConstitutiveModels.UniaxialStrain","text":"Provides an analytic motion for uniaxial strain\n\nThis is\n\nmathbfF = beginbmatrix\nlambda  0  0 \n0        1  0 \n0        0  1\nendbmatrix\n\n\n\n\n\n","category":"type"},{"location":"#ConstitutiveModels.cauchy_stress-Tuple{ConstitutiveModel, Vararg{Any, 5}}","page":"Home","title":"ConstitutiveModels.cauchy_stress","text":"non-AD Cauchy stress\n\n\n\n\n\n","category":"method"},{"location":"#ConstitutiveModels.deformation_gradient","page":"Home","title":"ConstitutiveModels.deformation_gradient","text":"Returns the deformation gradient for a given motion\n\n\n\n\n\n","category":"function"},{"location":"#ConstitutiveModels.deformation_gradient-Union{Tuple{T}, Tuple{Type{IsochoricUniaxialStress}, T}} where T<:Number","page":"Home","title":"ConstitutiveModels.deformation_gradient","text":"\n\n\n\n","category":"method"},{"location":"#ConstitutiveModels.deformation_gradient-Union{Tuple{T}, Tuple{Type{SimpleShear}, T}} where T<:Number","page":"Home","title":"ConstitutiveModels.deformation_gradient","text":"\n\n\n\n","category":"method"},{"location":"#ConstitutiveModels.deformation_gradient-Union{Tuple{T}, Tuple{Type{UniaxialStrain}, T}} where T<:Number","page":"Home","title":"ConstitutiveModels.deformation_gradient","text":"\n\n\n\n","category":"method"}]
}