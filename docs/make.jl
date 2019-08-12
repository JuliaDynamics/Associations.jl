using CausalityTools
using TimeseriesSurrogates
using PyCall, Conda
using Documenter, DocumenterMarkdown
#Conda.add("scipy")
using Plots
using DynamicalSystems
using Distributions
using StaticArrays
using Statistics
using StatsBase

ENV["GKSwstype"] = "100"

PAGES = [
    "index.md",
    "Syntax overview" => "syntax_overview.md",
    "CausalityToolsBase" => [
        "Discretization" => "causalitytoolsbase/discretization.md"
    ],
    "Transfer operator estimation" => [
        "Transfer operator" => "perronfrobenius/transferoperator.md",
        "Invariant measure" => "perronfrobenius/invariantmeasure.md"
    ],
    "Transfer entropy" => [
        "Estimators" => "transferentropy/transferentropy_estimators.md",
        "Convenience functions" => "transferentropy/convenience_functions_te.md"
    ]
]

makedocs(
    sitename = "CausalityTools.jl documentation",
    modules = [CausalityTools, TransferEntropy, PerronFrobenius, CrossMappings, CausalityToolsBase],
    format = DocumenterMarkdown.Markdown(),
    linkcheck = true,
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
