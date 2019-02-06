using CausalityTools
using TransferEntropy, PerronFrobenius, StateSpaceReconstruction, Simplices
using CrossMappings
using TimeseriesSurrogates
using Documenter, DocumenterMarkdown
using PyCall, Conda
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
    "Glossary" => [
        "Delay embeddings (SSR)" => "glossary/embed.md",
        "Discretization" => "glossary/discretization.md",
        "Invariantizing" => [
            "glossary/invariantizing/invariantizing.md",
            "glossary/invariantizing/invariantize.md",
            "glossary/invariantizing/forwardlinearmap_invariant.md",
        ]

    ],
    #"StateSpaceReconstruction" => [
        #"Delay embeddings" => "embed.md",
        #"Rectangular partitioning" => "partitioning_rectangular.md",
        #"Triangulation partitioning" => "partitioning_triangulation.md"
    #],

    "Plot recipes" => [
        "Recipe overview" => "plot_recipes/overview_recipes.md",
        "Transfer operator" => [
            "plot_recipes/recipe_transferoperator_triang_exact.md",
            "plot_recipes/recipe_transferoperator_triang_approx.md",
            "plot_recipes/recipe_transferoperator_rectangular_binning.md"
        ],
        "Invariant distribution" => [
            "plot_recipes/recipe_invariant_distribution.md"
        ]
    ],

    "PerronFrobenius" => [
        "Transfer operator" => [
			"transferoperator/overview.md",
            "transferoperator/transferoperator_rectangular_binning.md",
            "transferoperator/transferoperator_triang_exact.md",
            "transferoperator/transferoperator_triang_approx.md",
            "transferoperator/composite_types/triang_approx.md",
            "transferoperator/composite_types/triang_exact.md",
            "transferoperator/composite_types/rectangular_binning.md"
        #    "Transfer operator triangulation docs" => "transferoperator/transferoperator_triang_approx.md"
        ],
        #"Examples" => [
        #    "Transfer operator grid example" => "transferoperator/transferoperator_grid_example.md",
        #    "Transfer operator triangulation" => "transferoperator/transferoperator_triang_example.md"
        #],

        "Invariant measure" => [
			"invariantmeasure/invariantmeasure_overview.md",
            "invariantmeasure/invariantmeasure.md",
            "invariantmeasure/rectangularinvariantmeasure.md",
            "invariantmeasure/inducedinvariantmeasure.md",
            "invariantmeasure/averageinvariantmeasure.md"
        ],
    ],

    "Causality algorithms" => [
        "TransferEntropy" => [
            #"Common interface" => "transferentropy/commoninterface_TE.md",
            "kNN estimator" => "transferentropy/estimator_TE_kNN.md",
            "Transfer operator grid estimator" => "transferentropy/estimator_TE_transferoperator_grid.md",
            "Visitation frequency estimator" => "transferentropy/estimator_TE_visitfreq.md",
            "Wrapper" => "transferentropy/wrapper_TE.md"
        ],
        "Cross mappings" => [
            "crossmappings/ccm/overview.md",
            "crossmappings/ccm/crossmapping.md",
            "crossmappings/ccm/convergentcrossmapping.md"
        ]
    ],

    "Surrogate data" => [
        "The method of surrogate data" => "surrogates/surrogates_overview.md",
        "Docs" => [
            "randomshuffle" => "surrogates/randomshuffle_docs.md",
            "randomamplitudes" => "surrogates/randomamplitudes_docs.md",
            "randomphases" => "surrogates/randomphases_docs.md",
            "aaft" => "surrogates/aaft_docs.md",
            "iaaft" => "surrogates/iaaft_docs.md"
        ],
        "Examples" => [
            "randomshuffle" => "surrogates/randomshuffle_example.md",
            "randomamplitudes" => "surrogates/randomamplitudes_example.md",
            "randomphases" => "surrogates/randomphases_example.md",
            "aaft" => "surrogates/aaft_example.md",
            "iaaft" => "surrogates/iaaft_example.md"
        ]
    ],

    #"Workflow" => "workflow.md",

    #"Interop with DynamicalSystems.jl" => [
    #    "Example 1" => "interop_dynamicalsystems_infoflow1.md"
    #],

    "Examples" => [
        # Transfer operators
        "examples/transferoperators/exact.md",
        "examples/transferoperators/approx.md",

        # Cross mappings
        "examples/crossmappings/ccm_gif.md",
        "examples/crossmappings/examples_crossmappings_ar1.md",
        "examples/crossmappings/examples_crossmappings_henon2.md",
        "examples/crossmappings/examples_crossmappings_linearmap3d_nonlinearcoupling.md",
        "examples/crossmappings/examples_crossmappings_logistic2.md",
        "examples/crossmappings/examples_crossmappings_logistic3.md",
        "examples/crossmappings/examples_crossmappings_verdes.md"
    ],

    "Example systems" => [
		"ar1" => "example_systems/ar1.md",
        "henon2" => "example_systems/henon2.md",
        "logistic2" => "example_systems/logistic2.md",
        "logistic3" => "example_systems/logistic3.md",
		"linearmap3d_nonlinearcoupling" => "example_systems/linearmap3d_nonlinearcoupling.md",
		"verdes" => "example_systems/verdes.md"
    ]
]

makedocs(
    sitename = "CausalityTools.jl documentation",
    modules = [CausalityTools,
                TransferEntropy,
                PerronFrobenius,
                StateSpaceReconstruction,
                TimeseriesSurrogates,
                CrossMappings],
    format = Markdown(),
    pages = PAGES
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
