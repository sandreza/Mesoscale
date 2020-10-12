using Documenter

channel_simulations = [
    "Home" => "channels_home.md",
    "Simulation 1" => "simulation_1.md",
    "Abernathy Setup" => "simulation_abernathy.md",
]
makedocs(
    pages = [
        "Home" => "index.md",
        "Channels" => channel_simulations ,
    ],
    sitename = "Mesoscale",
    format = Documenter.HTML(collapselevel = 1),
)

deploydocs(
    repo = "github.com/sandreza/Mesoscale.git",
)