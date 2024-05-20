using Documenter, SmoothInterpolation

makedocs(;
    sitename = "SmoothInterpolation.jl",
    pages = [
        "Home" => "index.md",
        "Derivation" => [
            "Motivation " => "motivation.md",
            "Derivation of smoothed linear interpolation" => "derivation_smoothed_linear_interpolation.md",
        ],
    ],
    format = Documenter.HTML(;
        mathengine = Documenter.MathJax(
            Dict(:TeX => Dict(:equationNumbers => Dict(:autoNumber => "AMS"))),
        ),
    ),
)