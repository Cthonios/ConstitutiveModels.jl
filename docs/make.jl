using ConstitutiveModels
using Documenter

DocMeta.setdocmeta!(ConstitutiveModels, :DocTestSetup, :(using ConstitutiveModels); recursive=true)

makedocs(;
    modules=[ConstitutiveModels],
    authors="Craig M. Hamel <cmhamel32@gmail.com> and contributors",
    repo="https://github.com/cmhamel/ConstitutiveModels.jl/blob/{commit}{path}#{line}",
    sitename="ConstitutiveModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cmhamel.github.io/ConstitutiveModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/cmhamel/ConstitutiveModels.jl",
    devbranch="main",
)
