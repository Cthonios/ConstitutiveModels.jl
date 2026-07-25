using ConstitutiveModels
using Documenter
using LaTeXStrings
using RecipesBase
using Plots

DocMeta.setdocmeta!(ConstitutiveModels, :DocTestSetup, :(using ConstitutiveModels); recursive=true)

recipes_ext = Base.get_extension(ConstitutiveModels, :ConstitutiveModelsRecipesBaseExt)

makedocs(;
    modules=[ConstitutiveModels, recipes_ext],
    authors="Craig M. Hamel <cmhamel32@gmail.com> and contributors",
    repo="https://github.com/cmhamel/ConstitutiveModels.jl/blob/{commit}{path}#{line}",
    sitename="ConstitutiveModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cthonios.github.io/ConstitutiveModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home"               => "index.md",
        "Notation"           => "notation.md",
        "Abstract Interface" => "abstract_interface.md",
        "Models"             => [
            "models/finite_def_j2_plasticity.md"  
        ],
        "Modules"            => [
            "Heat conduction" => [
                "modules/heat_conduction/fouriers_law.md"
            ],
            "Hyperelasticity" => [
                "modules/hyperelasticity/hyperelasticity.md",
                "modules/hyperelasticity/gent.md",
                "modules/hyperelasticity/hencky.md",
                "modules/hyperelasticity/linear_elastic.md",
                "modules/hyperelasticity/mooney_rivlin.md",
                "modules/hyperelasticity/neohookean.md",
                "modules/hyperelasticity/saint_venant_kirchoff.md",
                "modules/hyperelasticity/seth_hill.md"
            ]  
        ],
        "Utils"               => [
            "Elastic Constants" => "utils/elastic_constants.md",
            "Material Symmetry" => "utils/material_symmetry.md",
            "Properties"        => "utils/properties.md",
            "Simple Motions"    => "utils/simple_motions.md"
        ]
        # "Mechanicals Models" => [
        #     "mechanics/mechanical_models.md",
        #     "mechanics/finite_def_j2_plasticity.md",
        #     "mechanics/gent.md",
        #     "mechanics/hencky.md",
        #     "mechanics/linear_elastic.md",
        #     "mechanics/mooney_rivlin.md",
        #     "mechanics/neohookean.md",
        #     "mechanics/saint_venant_kirchoff.md",
        #     "mechanics/seth_hill.md"
        # ],
        # "Thermal Models"     => "thermal_models.md",
        # "Common Methods"     => "common_methods.md",
        # "Elastic Constants"  => "elastic_constants.md",
        # "Simple Motions"     => "simple_motions.md"
    ],
    checkdocs=:none
)

deploydocs(;
    repo="github.com/Cthonios/ConstitutiveModels.jl",
    devbranch="main",
)
