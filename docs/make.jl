using Documenter, CachedInterpolations

makedocs(;
    sitename = "CachedInterpolations.jl",
    pages = [
        "Home" => "index.md",
        "Mathematical Construction" => [
            "Motivation " => "construction/motivation.md",
            "Construction of smoothed linear interpolation" => "construction/construction_smoothed_linear_interpolation.md",
            "Integrating" => "construction/integrating.md",
            "Inverting the integral" => "construction/inverting_the_integral.md",
        ],
        "Examples" => "examples.md",
        "To Cache or not to Cache" => "cache.md",
        "API" => "api.md",
    ],
    format = Documenter.HTML(;
        mathengine = Documenter.MathJax(
            Dict(:TeX => Dict(:equationNumbers => Dict(:autoNumber => "AMS"))),
        ),
    ),
)

deploydocs(; repo = "github.com/SouthEndMusic/CachedInterpolations.jl.git")
