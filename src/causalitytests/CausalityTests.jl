using Reexport

@reexport module CausalityTests
    using LaTeXStrings, RecipesBase
    import CausalityToolsBase: causality
    import NearestNeighbors
    import StatsBase
    import Distances
    import Distances: Chebyshev
    import CausalityToolsBase
    import CausalityToolsBase: CausalityTest, RectangularBinning
    import CrossMappings: crossmap, convergentcrossmap
    import TimeseriesSurrogates: randomshuffle
    import TransferEntropy: transferentropy, 
        TransferOperatorGrid, VisitationFrequency, NearestNeighbourMI, 
        BinningTransferEntropyEstimator
    import UncertainData
    import UncertainData:
        AbstractUncertainValue,
        AbstractUncertainValueDataset,
        AbstractUncertainIndexValueDataset,
        resample,
        ConstrainedResampling
    resample(v::Vector{Real}) = v

    import ..transferentropy

    include("causality.jl")

    ################################################################
    # Test definitions
    ################################################################
    # These are the basic tests
    include("tests/test_definitions.jl")

    # These meta tests apply the causality tests defined in 
    # `test_definitions.jl` in some more complicated manner, 
    # for example on random subsets of the data.
    include("MetaCausalityTest/MetaCausalityTest.jl")

    ################################################################
    # Test results
    ################################################################
    include("CausalAnalysis/causal_analyses.jl")

    ################################################################
    # On uncertain data with uncertainties in both index and value
    ################################################################
    include("causality_on_uncertaindata.jl")

    export causality
end



"""
	CausalityTests

A module defining causality tests.
"""
CausalityTests
