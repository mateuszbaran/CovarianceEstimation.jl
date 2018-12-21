using Documenter, CovarianceEstimation

makedocs(
    modules = [CovarianceEstimation],
    format = Documenter.HTML(),
    sitename = "CovarianceEstimation.jl",
    authors = "Mateusz Baran, Thibaut Lienart, and contributors.",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Methods" => "man/methods.md",
            "Linear shrinkage estimators" => "man/lshrink.md",
            "Nonlinear shrinkage estimators" => "man/nlshrink.md"
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md"
        ],
    ],
)
