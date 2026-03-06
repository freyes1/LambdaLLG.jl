using Pkg; Pkg.activate("../")

using Documenter
using LambdaLLG

makedocs(
    sitename = "LambdaLLG.jl",
    modules = [LambdaLLG],
    checkdocs = :exports,
    format = Documenter.HTML(prettyurls = get(ENV, "CI", "false") == "true"),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
)
