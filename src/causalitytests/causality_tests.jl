using Reexport

@reexport module CausalityTests
    import CausalityToolsBase: causality
    import NearestNeighbors
    import StatsBase
    import Distances
    import CausalityToolsBase
    import CausalityToolsBase: CausalityTest, RectangularBinning
    import CrossMappings: crossmap, convergentcrossmap
    import TimeseriesSurrogates: randomshuffle
    import TransferEntropy: transferentropy, TransferOperatorGrid, VisitationFrequency
    import UncertainData
    import UncertainData:
        AbstractUncertainValue,
        AbstractUncertainValueDataset,
        AbstractUncertainIndexValueDataset,
        resample,
        ConstrainedResampling


    resample(v::Vector{Real}) = v

    """
        causality(source::AbstractVector, target::AbstractVector, test::CausalityTest)

    Test for a causal influence from `source` to `target` using the provided causality `test`.
        
    ## Examples

    ```julia
    x, y = rand(300), rand(300)

    # Define some causality tests and apply them to `x` and `y`.
    test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 20)
    test_cm = CrossMappingTest(n_reps = 10)
    test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)
    test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1:5)
    test_jdd = JointDistanceDistributionTest()
    test_jddt = JointDistanceDistributionTTest()
    predtest = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = -5:5)
    test_pa = PredictiveAsymmetryTest(predictive_test = predtest)

    causality(x, y, test_ccm)
    causality(x, y, test_cm)
    causality(x, y, test_vf)
    causality(x, y, test_tog)
    causality(x, y, test_jdd)
    causality(x, y, test_jddt)
    causality(x, y, test_pa)
    ```
    """
    function causality(x::AbstractVector, y::AbstractVector, test::CausalityTest) end
    
    function summarise(test::CausalityTest)
        _type = typeof(test)

        strs = ["$fn = $(test.:($fn))" for fn in fieldnames(_type)]
        return "$_type" * "(" * join(strs, ", ") * ")"
    end

    Base.show(io::IO, test::CausalityTest) = print(io, summarise(test))


    ################################################################
    # Distance based causality tests
    ################################################################
    include("tests_distance_based/DistanceBasedCausalityTest.jl")

    # S-measure test
    # ---------------------
    import ..s_measure
    
    include("tests_distance_based/SMeasureTest.jl")
    
    # Joint distances tests
    # ---------------------
    import ..joint_distance_distribution

    include("tests_distance_based/JointDistancesCausalityTest.jl")
    include("tests_distance_based/JointDistanceDistributionTest.jl")
    include("tests_distance_based/JointDistanceDistributionTTest.jl")

    # Cross mapping tests
    # ---------------------
    include("tests_distance_based/CrossMappingTest.jl")
    include("tests_distance_based/ConvergentCrossMappingTest.jl")


    ################################################################
    # Entropy based causality tests
    ################################################################
    include("tests_entropy_based/EntropyBasedCausalityTest.jl")

    # Transfer entropy causality tests
    # ---------------------------------------
    include("tests_entropy_based/TransferEntropyCausalityTest.jl")
    include("tests_entropy_based/VisitationFrequencyTest.jl")
    include("tests_entropy_based/TransferOperatorGridTest.jl")
    include("tests_entropy_based/ApproximateSimplexIntersectionTest.jl")
    include("tests_entropy_based/ExactSimplexIntersectionTest.jl")

    ################################################################
    # Predictive asymmetry causality tests
    ################################################################
    include("tests_predictive_asymmetry/PredictiveAsymmetryTest.jl")

    ################################################################
    # On uncertain data with uncertainties in both index and value
    ################################################################
    include("causality_on_uncertaindata.jl")

    ################################################################
    # Preprocessed data + causality test
    ################################################################
    include("highlevel_tests/PreprocessedDataCausalityTest.jl")
    include("highlevel_tests/BinnedDataCausalityTest.jl")
    include("highlevel_tests/ConstrainedTest.jl")
    include("highlevel_tests/InterpolateBinTest.jl")
    include("highlevel_tests/RandomSequencesTest.jl")

    include("highlevel_tests/causality_BinnedDataCausalityTest.jl")
    include("highlevel_tests/causality_ConstrainedTest.jl")
    include("highlevel_tests/causality_InterpolateBinTest.jl")
    include("highlevel_tests/causality_RandomSequencesTest.jl")

    export causality
end



"""
	CausalityTests

A module defining causality tests.
"""
CausalityTests
