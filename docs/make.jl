using Documenter, ThermoCycleGlides

makedocs(
sitename="ThermoCycleGlides.jl",
modules = [ThermoCycleGlides],
doctest = true,
pages = [
    "Home" => "index.md",
    "Examples" => "examples.md",
    #"Optimization" => "optimization.md",
    # "Cycle Optimization" => "Optimization.md",
    "References" => "reference.md"
]
)

deploydocs(repo = "github.com/Sush1090/ThermoCycleGlides.jl",push_preview = true)