cd(@__DIR__)
using Pkg
CI = get(ENV, "CI", nothing) == "true" || get(ENV, "GITHUB_TOKEN", nothing) !== nothing
CI && Pkg.activate(@__DIR__)
CI && Pkg.instantiate()
ENV["GKSwstype"] = "100" # allow local builds without output
using Documenter
using DocumenterTools: Themes
using CausalityTools
using ComplexityMeasures
using StateSpaceSets

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
    "Correlation measures" => "correlation_measures.md",
    "Information measures" => [
        "Probability mass functions" => "probabilities.md",
        "Entropy" => "entropy.md",
        "Conditional entropy" => "entropy_conditional.md",
        "Mutual information" => "mutualinfo.md",
        "Conditional mutual information" => "condmutualinfo.md",
        "Transfer entropy" => "transferentropy.md",
    ],
    "Cross mappings" => "cross_mappings.md",
    "Joint distance distribution" => "jdd.md",
    "Independence testing" => "independence.md",
    "Examples" => [
        "Quickstart" => [
            "Mutual information" => "quickstart/quickstart_mi.md",
            "Independence testing" => "quickstart/quickstart_independence.md",
        ],
        "Longer examples" => [
            "Entropy" => "examples/examples_entropy.md",
            "Conditional entropy" => "examples/examples_conditional_entropy.md",
            "Mutual information" => "examples/examples_mutualinfo.md",
            "Transfer entropy" => "examples/examples_transferentropy.md",
            "Cross mappings" => "examples/examples_cross_mappings.md",
            "Independence testing" => "examples/examples_independence.md",
        ],
    ],
    "Experimental" => "experimental.md",
]

makedocs(
    modules = [CausalityTools, ComplexityMeasures, StateSpaceSets],
    format = Documenter.HTML(
        prettyurls = CI,
        sidebar_sitename = false,
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
