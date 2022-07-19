using ChebPade
using Documenter

DocMeta.setdocmeta!(ChebPade, :DocTestSetup, :(using ChebPade); recursive=true)

makedocs(;
    modules=[RobustChebPade],
    authors="Matthew Priddin and contributors",
    repo="https://github.com/mjp98/ChebPade.jl/blob/{commit}{path}#{line}",
    sitename="ChebPade.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mjp98.github.io/ChebPade.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mjp98/ChebPade.jl",
    devbranch="main",
)
