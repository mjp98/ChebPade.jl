using RobustChebPade
using Documenter

DocMeta.setdocmeta!(RobustChebPade, :DocTestSetup, :(using RobustChebPade); recursive=true)

makedocs(;
    modules=[RobustChebPade],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/RobustChebPade.jl/blob/{commit}{path}#{line}",
    sitename="RobustChebPade.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/RobustChebPade.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/RobustChebPade.jl",
    devbranch="main",
)
