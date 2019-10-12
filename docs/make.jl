using Documenter, UMAT

makedocs(;
    modules=[UMAT],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/JuliaFEM/UMAT.jl/blob/{commit}{path}#L{line}",
    sitename="UMAT.jl",
    authors="Tero Frondelius, Ivan Yashchuk",
    assets=String[],
)

deploydocs(;
    repo="github.com/JuliaFEM/UMAT.jl",
)
