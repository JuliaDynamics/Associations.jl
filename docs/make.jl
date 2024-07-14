cd(@__DIR__)
# Doc-specific (remaining packages are imported  in `build_docs_with_style.jl`, which is
# downloaded)
using DocumenterCitations
import Downloads

# Packages used in the doc build.
using CausalityTools
using ComplexityMeasures
using StateSpaceSets

pages = [
    "Overview" => "index.md",
    "Association measures" => "associations.md",
    "Independence" => "independence.md",
    "Network/graph inference" => "causal_graphs.md",
    "Examples" => [
        "Association measures" => "examples/examples_associations.md",
        #"Correlation examples" => "examples/examples_correlation.md",
        #"Closeness examples" => "examples/examples_closeness.md",
        #"Crossmap examples" => "examples/examples_cross_mappings.md",
    ],
    "Extended examples" => [
        "extended_examples/cross_mapping.md",
        "extended_examples/pairwise_asymmetric_inference.md",
        "extended_examples/mutual_information.md",
    ],
    "Basics and tutorials" => [
        "Encoding elements" => "encoding_tutorial.md",
        "Encoding input datasets" => "discretization_tutorial.md",
        "Counts and probabilities" => "probabilities_tutorial.md",
        "Information measures" => "info_tutorial.md",

    ],
    # "Independence testing" => "independence.md",
    # "Causal graphs" => "causal_graphs.md",
    # "Predefined systems" => "coupled_systems.md",
    # "Experimental" => "experimental.md",
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
