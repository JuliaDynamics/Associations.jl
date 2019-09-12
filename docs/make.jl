using CausalityTools
using TimeseriesSurrogates
using PyCall, Conda
using HypothesisTests
using Documenter, DocumenterMarkdown
#Conda.add("scipy")
using Plots
using DynamicalSystems
using Distributions
using StaticArrays
using Statistics
using StatsBase
using UncertainData

ENV["GKSwstype"] = "100"

PAGES = [
    "index.md",
    "Syntax overview" => "syntax_overview.md",
    "CausalityToolsBase" => [
        "Discretization" => "causalitytoolsbase/discretization.md",
        "Delay reconstructions" => "causalitytoolsbase/delay_reconstructions.md"
    ],
    "Transfer operator estimation" => [
        "Transfer operator" => "perronfrobenius/transferoperator.md",
        "Invariant measure" => "perronfrobenius/invariantmeasure.md"
    ],
    "Transfer entropy" => [
        "Estimators" => "transferentropy/transferentropy_estimators.md",
        "Convenience functions" => "transferentropy/convenience_functions_te.md",
        "TEVars" => "transferentropy/TEVars.md"
    ],
    "Distance based measures" => [
        "CCM" => [
            "Overview" => "crossmappings/ccm/overview.md",
            "Cross mapping" => "crossmappings/ccm/crossmapping.md",
            "Converent cross mapping" => "crossmappings/ccm/convergentcrossmapping.md"
        ], 
        "Joint distance distribution" => "algorithms/joint_distance_distribution.md"
    ],
    "Examples" => [
        "CCM" => [
            "Henon map" => "examples/crossmappings/examples_crossmappings_henon2.md"
        ]
    ]
]

makedocs(
    sitename = "CausalityTools.jl documentation",
    modules = [CausalityTools, TransferEntropy, PerronFrobenius, CrossMappings, CausalityToolsBase, UncertainData],
    format = DocumenterMarkdown.Markdown(),
    linkcheck = false,
    pages = PAGES,
    highlightsig = true
)

if !Sys.iswindows()
    deploydocs(
        deps   = Deps.pip("mkdocs==0.17.5", "mkdocs-material==2.9.4",
        "python-markdown-math", "pygments", "pymdown-extensions"),
        repo   = "github.com/kahaaga/CausalityTools.jl.git",
        target = "site",
        make = () -> run(`mkdocs build`)
    )
end
