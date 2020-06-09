using Pkg
# cd(@__DIR__)
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
# CI && Pkg.activate(@__DIR__)
# CI && Pkg.instantiate()
# CI && (ENV["GKSwstype"] = "100")
using CausalityToolsBase
using PerronFrobenius
using TransferEntropy
using CrossMappings
using CausalityTools
using DynamicalSystems
using Plots
using LaTeXStrings
using StaticArrays

using Simplices
using Documenter
using DocumenterTools: Themes

# %% JuliaDynamics theme.
# download the themes
using DocumenterTools: Themes
for file in ("juliadynamics-lightdefs.scss", "juliadynamics-darkdefs.scss", "juliadynamics-style.scss")
    download("https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/$file", joinpath(@__DIR__, file))
end
# create the themes
for w in ("light", "dark")
    header = read(joinpath(@__DIR__, "juliadynamics-style.scss"), String)
    theme = read(joinpath(@__DIR__, "juliadynamics-$(w)defs.scss"), String)
    write(joinpath(@__DIR__, "juliadynamics-$(w).scss"), header*"\n"*theme)
end
# compile the themes
Themes.compile(joinpath(@__DIR__, "juliadynamics-light.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-light.css"))
Themes.compile(joinpath(@__DIR__, "juliadynamics-dark.scss"), joinpath(@__DIR__, "src/assets/themes/documenter-dark.css"))

# %% Build docs
#cd(@__DIR__)
#ENV["JULIA_DEBUG"] = "Documenter"

PAGES = [
    "Overview" => "index.md",
    "Transfer entropy" => "TransferEntropy.md",
    "PredictiveAsymmetry" => "PredictiveAsymmetry.md",
    "Joint distance distribution" => "JointDistanceDistribution.md",
    "S-measure" => "SMeasure.md",
    "Cross mapping" => "CrossMapping.md",
    "Example applications" => [
        "Transfer entropy" => "example_applications/example_transferentropy.md",
        "Predictive asymmetry" => "example_applications/example_predictive_asymmetry.md"
    ],
    "Pre-defined coupled systems" => "CoupledSystems.md",
    "Utility tools" => [
        "CausalityToolsBase" => "CausalityToolsBase.md",
        "PerronFrobenius" => "PerronFrobenius.md"
    ]
]

makedocs(
    modules = [CausalityToolsBase, PerronFrobenius, TransferEntropy, CrossMappings, CausalityTools],
    format = Documenter.HTML(
        prettyurls = CI,
        assets = [
            asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
        ),
    sitename = "CausalityTools.jl",
    authors = "Kristian Agasøster Haaga, David Diego, Tor Einar Møller, Bjarte Hannisdal, George Datseris",
    pages = PAGES
)

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/CausalityTools.jl.git",
        target = "build",
        push_preview = true
    )
end
