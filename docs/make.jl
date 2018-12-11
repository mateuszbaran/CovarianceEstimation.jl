using Documenter, CovarianceEstimation

makedocs(
    modules = [CovarianceEstimation],
    format = :html,
    html_prettyurls = false,
    sitename = "CovarianceEstimation.jl",
    authors = "Mateusz Baran, Thibaut Lienart, and contributors.",
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Methods" => "man/methods.md"
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md"
        ],
    ],
)
