using ConstitutiveModels

inputs = Dict(
    "FeFv" => Dict(
        "equilibrium branch" => Dict(
            "strain energy density" => "NeoHookean",
            "bulk modulus" => 1000.0,
            "shear modulus" => 1.0
        ),
        "maxwell branches" => [
            Dict(
                "shear modulus" => 10.0,
                "relaxation time" => 1.0
            )
        ]
    )
)

model, props, Z_old, Z_new = ConstitutiveModels.initialize(inputs, "FeFv")
