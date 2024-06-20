using Feynman
using Documenter

DocMeta.setdocmeta!(Feynman, :DocTestSetup, :(using Feynman); recursive=true)

makedocs(;
    modules=[Feynman],
    clean = true,
    checkdocs = :none,
    doctest = true,
    authors = "Dushan Priyasad, Janko Boehm",
    sitename = "Feynman.jl",
    expandfirst = ["Overview.md"],
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://singular-gpispace.github.io/Feynman.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => [
            "index.md",
            "Installation.md",
         #   "SmallExample.md",
        ],

        "Examples" =>[
            "FeynmanIBP.md",

        ],

        "Functions" => "Overview.md",

    ],


)

deploydocs(;
    repo="github.com/singular-gpispace/Feynman.jl",
    devbranch="main",
)
