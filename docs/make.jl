push!(LOAD_PATH, "/Users/antoinerestivo/Desktop/Permanents.jl-main/src")

using Documenter, Permanents

# Functions = "Functions" => "index.md"

# About = "Introduction" => "index.md"
#
# GettingStarted = "gettingstarted.md"
#
# Examples = "Examples" => [
#         "examples/flux.md"
#     ]
#
# License = "License" => "license.md"

# PAGES = [Functions]

makedocs(
    sitename = "Permanents.jl",
    modules = [Permanents],
    authors = "Benoit Seron, Antoine Restivo",
    format = Documenter.HTML(),
)

# makedocs(
#     sitename = "Permanents.jl",
#     modules = [Permanents],
#     authors = "Benoit Seron, Antoine Restivo",
#     format = Documenter.HTML(),
#     pages = PAGES
# )

# deploydocs(repo = "github.com/Evizero/Augmentor.jl.git")
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
