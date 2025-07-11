module ConstitutiveModelsRecipesBaseExt

using ConstitutiveModels
using LaTeXStrings
using RecipesBase

@recipe function f(
    ::PureShearStrain, 
    ∇us, σs, Zs,
    σ_11s_an=nothing, 
    σ_33s_an=nothing, 
    σ_12s_an=nothing
)
    # check sizes here
    @assert length(∇us) == length(σs)
    @assert length(∇us) == length(Zs)
    if σ_11s_an !== nothing
        @assert length(∇us) == length(σ_11s_an)
    end
    if σ_33s_an !== nothing
        @assert length(∇us) == length(σ_33s_an)
    end
    if σ_12s_an !== nothing
        @assert length(∇us) == length(σ_12s_an)
    end

    # setup up global plot stuff
    xlabel --> "Shear"
    ylabel --> "Cauchy Stress"
    
    # unpack/re-format inputs
    Fs = map(x -> x + one(x), ∇us)
    F_12s = map(x -> x[1, 2], Fs)
    σ_11s = map(x -> x[1, 1], σs)
    σ_33s = map(x -> x[3, 3], σs)
    σ_12s = map(x -> x[1, 2], σs)

    num_pts = 10
    step = length(∇us) ÷ num_pts

    @series begin
        color      := :black
        label      := L"\mathbf{\sigma}_{11}"
        linewidth  := 3
        seriestype := :line
        F_12s, σ_11s
    end

    if σ_11s_an !== nothing
        @series begin
            color      := :black
            label      := L"\mathbf{\sigma}_{11} - Analytic"
            seriestype := :scatter
            F_12s[1:step:end], σ_11s_an[1:step:end]
        end
    end
    @series begin
        color      := :red
        label      := L"\mathbf{\sigma}_{33}"
        linewidth  := 3
        seriestype := :line
        F_12s, σ_33s
    end

    if σ_33s_an !== nothing
        @series begin
            color      := :red
            label      := L"\mathbf{\sigma}_{33} - Analytic"
            seriestype := :scatter
            F_12s[1:step:end], σ_33s_an[1:step:end]
        end
    end

    @series begin
        color      := :blue
        label      := L"\mathbf{\sigma}_{12}"
        linewidth  := 3
        seriestype := :line
        F_12s, σ_12s
    end

    if σ_12s_an !== nothing
        @series begin
            color      := :blue
            label      := L"\mathbf{\sigma}_{12} - Analytic"
            seriestype := :scatter
            F_12s[1:step:end], σ_12s_an[1:step:end]
        end
    end
end

@recipe function f(
    ::SimpleShear, 
    ∇us, σs, Zs,
    σ_11s_an=nothing, 
    σ_22s_an=nothing, 
    σ_12s_an=nothing
)
    # check sizes here
    @assert length(∇us) == length(σs)
    @assert length(∇us) == length(Zs)
    if σ_11s_an !== nothing
        @assert length(∇us) == length(σ_11s_an)
    end
    if σ_22s_an !== nothing
        @assert length(∇us) == length(σ_22s_an)
    end
    if σ_12s_an !== nothing
        @assert length(∇us) == length(σ_12s_an)
    end

    # setup up global plot stuff
    xlabel --> "Shear"
    ylabel --> "Cauchy Stress"
    
    # unpack/re-format inputs
    Fs = map(x -> x + one(x), ∇us)
    F_12s = map(x -> x[1, 2], Fs)
    σ_11s = map(x -> x[1, 1], σs)
    σ_22s = map(x -> x[2, 2], σs)
    σ_12s = map(x -> x[1, 2], σs)

    num_pts = 10
    step = length(∇us) ÷ num_pts

    @series begin
        color      := :black
        label      := L"\mathbf{\sigma}_{11}"
        linewidth  := 3
        seriestype := :line
        F_12s, σ_11s
    end

    if σ_11s_an !== nothing
        @series begin
            color      := :black
            label      := L"\mathbf{\sigma}_{11} - Analytic"
            seriestype := :scatter
            F_12s[1:step:end], σ_11s_an[1:step:end]
        end
    end

    @series begin
        color      := :red
        label      := L"\mathbf{\sigma}_{22}"
        linewidth  := 3
        seriestype := :line
        F_12s, σ_22s
    end

    if σ_22s_an !== nothing
        @series begin
            color      := :red
            label      := L"\mathbf{\sigma}_{22} - Analytic"
            seriestype := :scatter
            F_12s[1:step:end], σ_22s_an[1:step:end]
        end
    end

    @series begin
        color      := :blue
        label      := L"\mathbf{\sigma}_{12}"
        linewidth  := 3
        seriestype := :line
        F_12s, σ_12s
    end

    if σ_12s_an !== nothing
        @series begin
            color      := :blue
            label      := L"\mathbf{\sigma}_{12} - Analytic"
            seriestype := :scatter
            F_12s[1:step:end], σ_12s_an[1:step:end]
        end
    end
end

@recipe function f(
    ::M, 
    ∇us, σs, Zs,
    σ_11s_an=nothing, 
    σ_22s_an=nothing
) where M <: Union{
    <:UniaxialStrain,
    <:UniaxialStressDisplacementControl
}
    # check sizes here
    @assert length(∇us) == length(σs)
    @assert length(∇us) == length(Zs)
    if σ_11s_an !== nothing
        @assert length(∇us) == length(σ_11s_an)
    end
    if σ_22s_an !== nothing
        @assert length(∇us) == length(σ_22s_an)
    end

    # setup up global plot stuff
    xlabel --> "Stretch"
    ylabel --> "Cauchy Stress"
    
    # unpack/re-format inputs
    Fs = map(x -> x + one(x), ∇us)
    F_11s = map(x -> x[1, 1], Fs)
    σ_11s = map(x -> x[1, 1], σs)
    σ_22s = map(x -> x[2, 2], σs)

    num_pts = 10
    step = length(∇us) ÷ num_pts

    @series begin
        color      := :black
        label      := L"\mathbf{\sigma}_{11}"
        linewidth  := 3
        seriestype := :line
        F_11s, σ_11s
    end

    if σ_11s_an !== nothing
        @series begin
            color      := :black
            label      := L"\mathbf{\sigma}_{11} - Analytic"
            seriestype := :scatter
            F_11s[1:step:end], σ_11s_an[1:step:end]
        end
    end

    @series begin
        color      := :red
        label      := L"\mathbf{\sigma}_{22}"
        linewidth  := 3
        seriestype := :line
        F_11s, σ_22s
    end

    if σ_22s_an !== nothing
        @series begin
            color      := :red
            label      := L"\mathbf{\sigma}_{22} - Analytic"
            seriestype := :scatter
            F_11s[1:step:end], σ_22s_an[1:step:end]
        end
    end
end

end # module
