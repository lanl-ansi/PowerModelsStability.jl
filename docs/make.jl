using Documenter
using PowerModelsStability

makedocs(
    sitename = "PowerModelsStability",
    format = Documenter.HTML(),
    modules = [PowerModelsStability]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
