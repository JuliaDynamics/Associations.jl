cd(@__DIR__)
# Doc-specific (remaining packages are imported  in `build_docs_with_style.jl`, which is
# downloaded)
using DocumenterCitations
using DocumenterInterLinks
import Downloads

# Packages used in the doc build.
using CausalityTools
using ComplexityMeasures
using StateSpaceSets

pages = [
    "CausalityTools.jl" => "index.md",
    "Association measures" => "associations.md",
    "Independence" => "independence.md",
    "Network/graph inference" => "causal_graphs.md",
    "Examples" => [
        "Associations" => "examples/examples_associations.md", 
        "Independence testing" => "examples/examples_independence.md",
        "Causal graph inference" => "examples/examples_infer_graphs.md",
    ],
   "Tutorials" => [
        "Discretization" => "tutorials/discretization.md",
        "Counts and probabilities" => "tutorials/probabilities_tutorial.md",
        "Information measures" => "tutorials/info_tutorial.md",
     ],
    "References" => "references.md",
]


# Downloads.download(
#     "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
#     joinpath(@__DIR__, "build_docs_with_style.jl")
# )
include("build_docs_with_style.jl")

bibliography = CitationBibliography(
    joinpath(@__DIR__, "refs.bib");
    style=:authoryear
)

build_docs_with_style(pages, CausalityTools, ComplexityMeasures, StateSpaceSets;
    expandfirst = ["index.md"],
    bib = bibliography,
    pages = pages,
    authors = "Kristian Agasøster Haaga, David Diego, Tor Einar Møller, George Datseris",
)
