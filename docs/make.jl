using CausalityTools
using TimeseriesSurrogates
using PyCall, Conda
using HypothesisTests
using Documenter, DocumenterMarkdown
#Conda.add("scipy")
using Plots
using DynamicalSystems
using SimpleDiffEq
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
    "CHANGELOG.md",
    "Causality tests" => [
        "causalitytests/causality_from_time_series.md",
        "causalitytests/causality_from_dynamical_systems.md",
        "causalitytests/causality_tests.md",
        "causalitytests/CausalityTest.md",
        "causalitytests/ConvergentCrossMappingTest.md",
        "causalitytests/CrossMappingTest.md",
        "causalitytests/DistanceBasedCausalityTest.md",
        "causalitytests/EntropyBasedCausalityTest.md",
        "causalitytests/JointDistanceDistributionTest.md",
        "causalitytests/PredictiveAsymmetryTest.md",
        "causalitytests/TransferEntropyTest.md",
        "causalitytests/TransferOperatorGridTest.md",
        "causalitytests/VisitationFrequencyTest.md",
        "causalitytests/ApproximateSimplexIntersectionTest.md",
        "causalitytests/ExactSimplexIntersectionTest.md",

        "causalitytests/causality_from_uncertain_data.md",
        "causalitytests/BinnedDataCausalityTest.md"
    ],
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
        "TEVars" => "transferentropy/TEVars.md",
        "Effect of discretization scheme" => "transferentropy/examples_TE_different_partitionings.md",
    ],
    "Distance based measures" => [
        "CCM" => [
            "Overview" => "crossmappings/ccm/overview.md",
            "Cross mapping" => "crossmappings/ccm/crossmapping.md",
            "Converent cross mapping" => "crossmappings/ccm/convergentcrossmapping.md"
        ], 
        "Joint distance distribution" => "algorithms/joint_distance_distribution.md"
    ],
    "Worked example" => [
        "worked_examples/worked_example_transferentropy.md"
    ],
    
    "Example systems" => [
        "example_systems/noise.md",
        "example_systems/example_systems_overview.md",
        "example_systems/example_systems_discrete.md",
        "example_systems/example_systems_continuous.md"
    ],

    "Surrogate data" => [
        "surrogates/iaaft_docs.md",
        "surrogates/aaft_docs.md",
        "surrogates/randomphases_docs.md",
        "surrogates/randomamplitudes_docs.md",
        "surrogates/randomshuffle_docs.md",
        "surrogates/surrogates_overview.md"
    ],

    "Tutorials" => [
        "causality" => [
            "tutorials/list_of_tutorials.md",
            "tutorials/causality/binned_uncertain_data/tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedResampling.md",
            "tutorials/causality/binned_uncertain_data/tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedMeanResampling.md"
        ]
    ]
]

makedocs(
    sitename = "CausalityTools.jl documentation",
    #modules = [CausalityTools, TransferEntropy, PerronFrobenius, CrossMappings, CausalityToolsBase, UncertainData],
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
