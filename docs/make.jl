using Documenter

makedocs(
    sitename = "Mesoscale",
    format = Documenter.HTML(),
)

deploydocs(
    repo = "github.com/sandreza/Mesoscale.git",
)