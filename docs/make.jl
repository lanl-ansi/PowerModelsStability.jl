using Documenter
using PowerModelsStability

makedocs(
    sitename = "PowerModelsStability",
    authors = "Haoxiang Yang, Harsha Nagarajan, David M Fobes",
    format = Documenter.HTML(
        analytics = "",
        mathengine = Documenter.MathJax(),
        prettyurls = false,
        collapselevel = 1,
    ),
    modules = [PowerModelsStability],
    strict = false,
    pages = [
        "Introduction" => "index.md",
        "installation.md",
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/lanl-ansi/PowerModelsStability.jl.git",
    devbranch = "main",
    push_preview = true
)
