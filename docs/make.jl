using SPEX
using Documenter

DocMeta.setdocmeta!(SPEX, :DocTestSetup, :(using SPEX); recursive=true)

makedocs(;
    modules=[SPEX],
    authors="Will Kimmerer <kimmerer@mit.edu> and contributors",
    repo="https://github.com/"Wimmerer"/SPEX.jl/blob/{commit}{path}#{line}",
    sitename="SPEX.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://"Wimmerer".github.io/SPEX.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/"Wimmerer"/SPEX.jl",
    devbranch="main",
)
