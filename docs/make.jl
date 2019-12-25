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
    #"index.md",
    #"FAQ.md",
    "overview.md",
    "Syntax overview" => "syntax_overview.md",
    "CHANGELOG.md",
    "Causality tests" => [
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
        "causalitytests/abstract_test_types.md",
        "causalitytests/causality_tests_overview.md",
        "causalitytests/ConvergentCrossMappingTest.md",
        "causalitytests/CrossMappingTest.md",
        "causalitytests/JointDistanceDistributionTest.md",
        "causalitytests/SMeasureTest.md",

        "causalitytests/NearestNeighbourMITest.md",
        "causalitytests/TransferOperatorGridTest.md",
        "causalitytests/VisitationFrequencyTest.md",
        "causalitytests/ApproximateSimplexIntersectionTest.md",
        "causalitytests/ExactSimplexIntersectionTest.md",

        "causalitytests/PredictiveAsymmetryTest.md",
        "causalitytests/NormalisedPredictiveAsymmetryTest.md"
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
        "Overview" => "transferentropy/overview_te.md",
        "Estimators" => "transferentropy/TE_estimators.md",
        "Estimating TE" => "transferentropy/estimating_TE.md",
        "Assigning marginals" => "transferentropy/assigning_marginals.md",
        "Generalised delay embedding" => "transferentropy/generalised_delay_reconstructions.md",
        "Convenience methods" => "transferentropy/convenience_methods_te.md",
        "Effect of discretization scheme" => "transferentropy/examples_TE_different_partitionings.md",
    ],
    "PredictiveAsymmetry" => [
        "Predictive asymmetry" => "PredictiveAsymmetry/predictive_asymmetry.md"
    ],
    "Distance based statistics" => [
        "CCM" => [
            "Overview" => "crossmappings/ccm/overview.md",
            "Cross mapping" => "crossmappings/ccm/crossmapping.md",
            "Converent cross mapping" => "crossmappings/ccm/convergentcrossmapping.md"
        ]
    ],

    "Joint distance distribution" => [
        "JointDistanceDistribution/joint_distance_distribution.md"
    ],

    "SMeasure" => [
        "S-measure" => "SMeasure/s_measure.md"
    ],
    "Worked example" => [
        "worked_examples/worked_example_transferentropy.md"
    ],

    "Simplices" => [
        "simplices/Simplices.md"
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
            "tutorials/causality/Tutorial_BinnedResampling_RandomSequencesTest_PredictiveAsymmetryTest_ar1_unidir.md",
            "tutorials/causality/binned_uncertain_data/tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedResampling.md",
            "tutorials/causality/binned_uncertain_data/tutorial_BinnedDataCausalityTest_PredictiveAsymmetryTest_BinnedMeanResampling.md"
        ],
        "tutorials/causality/tutorial_SMeasureTest.md"
    ]
]


# Some things need to be expanded before index.md for headers to work properly.
expandfirst = [
    "causalitytoolsbase/discretization.md",
    "causalitytoolsbase/delay_reconstructions.md",
    "simplices/Simplices.md",
    "causalitytests/causality_tests_overview.md",
    "perronfrobenius/transferoperator.md",
    "perronfrobenius/invariantmeasure.md",
    "causalitytests/NearestNeighbourMITest.md"
]

makedocs(
    sitename = "CausalityTools.jl documentation",
    #modules = [CausalityTools, TransferEntropy, PerronFrobenius, CrossMappings, CausalityToolsBase, UncertainData],
    modules = [CausalityTools, TransferEntropy, PerronFrobenius, CrossMappings, CausalityToolsBase, UncertainData, Simplices],
    format = DocumenterMarkdown.Markdown(),
    linkcheck = false,
    clean = true,
    expandfirst = expandfirst,
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
