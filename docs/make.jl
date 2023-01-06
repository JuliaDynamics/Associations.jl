cd(@__DIR__)
using Pkg
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
CI && Pkg.activate(@__DIR__)
CI && Pkg.instantiate()
ENV["GKSwstype"] = "100" # allow local builds without output
using Documenter
using DocumenterTools: Themes
using ComplexityMeasures
using CausalityTools

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
cd(@__DIR__)
ENV["JULIA_DEBUG"] = "Documenter"

PAGES = [
    "Overview" => "index.md",

    "Information measures" => [
        "Design philosophy" => "information_basics.md",
        "Contingency tables" => "contingency.md",
        "Probabilities" => "probabilities.md",
        "Entropy" => "entropy.md",

        "Conditional entropy" => "entropy_conditional.md",
        "Mutual information" => "mutualinfo.md",

        "Conditional mutual information" => "condmutualinfo.md",
        "Examples" => [
            "Entropy" => "examples/examples_entropy.md",
            "Mutual information" => "examples/examples_mutualinfo.md"
        ],
        "Timeseries analysis" => [
            "Transfer entropy" => "transferentropy.md",
            # "predictive_asymmetry.md",
        ]
    ],
    "Independence tests" => "independence.md",

    "Cross mappings" => "cross_mappings.md",
    "Joint distance distribution" => "jdd.md",
]

makedocs(
    modules = [CausalityTools, ComplexityMeasures, StateSpaceSets],
    format = Documenter.HTML(
        prettyurls = CI,
        assets = [
            asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
        ),
    sitename = "CausalityTools.jl",
    authors = "Kristian Agasøster Haaga, David Diego, Tor Einar Møller, George Datseris",
    pages = PAGES
)

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/CausalityTools.jl.git",
        target = "build",
        push_preview = true
    )
end
