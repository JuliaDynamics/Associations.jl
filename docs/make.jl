cd(@__DIR__)
# Doc-specific (remaining packages are imported  in `build_docs_with_style.jl`, which is
# downloaded)
using DocumenterCitations
using DocumenterInterLinks
import Downloads

# Packages used in the doc build.
using Associations
using ComplexityMeasures
using StateSpaceSets

pages = [
    "Associations.jl" => "index.md",
    "Core API reference" => [
        "Association measures" => "associations.md",
        "Independence" => "independence.md",
        "Network/graph inference" => "causal_graphs.md",
    ],
    "Extended API reference" => [
        "api/discretization_counts_probs_api.md",
        "api/counts_and_probabilities_api.md",
        "api/information_single_variable_api.md",
        "api/information_multivariate_api.md",
        "api/cross_map_api.md",
    ],
    "Examples" => [
        "Associations" => "examples/examples_associations.md",
        "Independence testing" => "examples/examples_independence.md",
        "Causal graph inference" => "examples/examples_infer_graphs.md",
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

build_docs_with_style(pages, Associations, StateSpaceSets;
    expandfirst=["index.md"],
    bib=bibliography,
    pages=pages,
    authors="Kristian Agasøster Haaga, David Diego, Tor Einar Møller, George Datseris",
)
