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
        # "Mechanical Models"  => "mechanical_models.md",
        "Mechanicals Models" => [
            "mechanics/mechanical_models.md",
            "mechanics/gent.md",
            "mechanics/hencky.md",
            "mechanics/linear_elastic.md",
            "mechanics/neohookean.md"
        ],
        "Thermal Models"     => "thermal_models.md",
        "Abstract Interface" => "abstract_interface.md",
        "Common Methods"     => "common_methods.md",
        "Elastic Constants"  => "elastic_constants.md",
        "Simple Motions"     => "simple_motions.md"
    ],
)

deploydocs(;
    repo="github.com/Cthonios/ConstitutiveModels.jl",
    devbranch="main",
)
