using Documenter, CausalityTools
push!(LOAD_PATH,"../src/")
ENV["GKSwstype"] = "100"

PAGES = [
    "Overview" => "index.md",
	"State space reconstruction" => [
	#	"Tutorial" => "reconstruction/embed_tutorial.md",
	"Quickstart" => "reconstruction/embed_quickstart.md"
		#"Partitioning" => "reconstruction/partitioning.md"
	],
    "Transfer (Perron-Frobenius) operator" => "transferoperator/transferoperator.md",
    "Probability estimation and invariant measures" => [
        "Invariant measures" => "invariantmeasure/invariantmeasure.md"
    ],
    "Surrogate data" => [
        "What is a surrogate?" => "surrogates/surrogates.md",
        "Constrained surrogate examples" => [
            "Random shuffle surrogates" => "surrogates/random_shuffle.md",
            "AAFT" => "surrogates/aaft.md",
            "IAAFT" => "surrogates/iaaft.md"
        ],
        "Unconstrained surrogate examples" => [
            "Random phase Fourier surrogates" => "surrogates/Fourierphase.md",
            "Random amplitude Fourier surrogates" => "surrogates/Fourieramp.md"
        ],
        "Function documentation" => "surrogates/functiondocs.md"
    ],
    "Causality measures" => [
        "causalitytools.md",
        "Transfer entropy" => [
            "Which estimator to choose?" => "transferentropy/overview.md",
            "Transfer operator estimator" => "transferentropy/transferoperator.md",
            "Visitation frequency estimator" => "transferentropy/visitfreq.md",
            "Quick start examples" => "transferentropy/quickstart.md",
            "Reducing bias" => "transferentropy/reducing_bias.md"
            ]
        ],
    "Examples of coupled dynamical systems" => [
        "overview" => "examplesystems/examples.md",
        "Discrete maps" => "examplesystems/discrete_maps.md",
        "Continuous" => "examplesystems/continuous.md"
    ],
    "Tutorials" => "tutorials/tutorialoverview.md"
]

makedocs(
    modules = [CausalityTools],
    format = :markdown,
    sitename = "CausalityTools.jl",
    authors = "Kristian Agasøster Haaga",
    pages = PAGES,
    clean = true
    # Use clean URLs, unless built as a "local" build
    #html_prettyurls = !("local" in ARGS),
    #html_canonical = "https://kahaaga.github.io/CausalityTools.jl/latest/"
)

# deploydocs(
#     Deps.pip("pygments", "mkdocs", "python-markdown-math")
#     repo   = "github.com/kahaaga/CausalityTools.jl.git",
#     julia  = "0.6",
#     target = "build",
#     deps = nothing,
#     make = nothing,
#     osname = "linux"
# )
