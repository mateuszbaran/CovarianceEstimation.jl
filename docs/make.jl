using Documenter, CovarianceEstimation, StatsBase

makedocs(
    modules = [CovarianceEstimation],
    format = Documenter.HTML(; prettyurls=get(ENV, "CI", "false") == "true"),
    sitename = "CovarianceEstimation.jl",
    authors = "Mateusz Baran, Thibaut Lienart, and contributors.",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Methods" => "man/methods.md",
            "Linear shrinkage estimators" => "man/lshrink.md",
            "Nonlinear shrinkage estimators" => "man/nlshrink.md",
            "MSE comparisons" => "man/msecomp.md"
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md"
        ],
    ],
)

deploydocs(
    repo = "github.com/mateuszbaran/CovarianceEstimation.jl.git", push_preview = true
)
