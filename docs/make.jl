using CausalityTools
using CausalityToolsBase
using TimeseriesSurrogates
using PerronFrobenius
using PyCall, Conda
using CrossMappings
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
using TransferEntropy
using UncertainData

ENV["GKSwstype"] = "100"

PAGES = [
    "index.md",
    "Syntax overview" => "syntax_overview.md",
    "CHANGELOG.md",
    "Causality tests" => [
        # From time series
        "causalitytests/causality_from_time_series.md",

        # From uncertain data
        "causalitytests/causality_from_uncertain_data_naive.md",
        "causalitytests/causality_from_uncertain_data_binneddatacausalitytest.md",
        "causalitytests/causality_from_uncertain_data_strictlyincreasing_interpolated.md",
        "causalitytests/causality_from_uncertain_data_naive_constrained.md",
        "causalitytests/causality_from_uncertain_data_InterpolateBinTest.md",
        "causalitytests/causality_from_uncertain_data_RandomSequencesTest.md",
        "causalitytests/causality_from_uncertain_data.md",

        # From dynamical systems
        "causalitytests/causality_from_dynamical_systems.md",

        # Tests
        "causalitytests/causality_tests.md",

        "causalitytests/ConvergentCrossMappingTest.md",
        "causalitytests/CrossMappingTest.md",
        "causalitytests/JointDistanceDistributionTest.md",
        "causalitytests/SMeasureTest.md",

        "causalitytests/TransferOperatorGridTest.md",
        "causalitytests/VisitationFrequencyTest.md",
        "causalitytests/ApproximateSimplexIntersectionTest.md",
        "causalitytests/ExactSimplexIntersectionTest.md",

        "causalitytests/PredictiveAsymmetryTest.md"
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
    "UncertainData" => [
        "uncertaindata/BinnedDataCausalityTest.md"
    ],
    "Tutorials" => [
        "causality" => [
            "tutorials/list_of_tutorials.md",
            "tutorials/causality/binned_uncertain_data/tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedResampling.md",
            "tutorials/causality/binned_uncertain_data/tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedMeanResampling.md"
        ],
        "tutorials/tutorial_SMeasureTest.md"
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
