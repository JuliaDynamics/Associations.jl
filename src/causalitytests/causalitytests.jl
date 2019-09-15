using Reexport

@reexport module CausalityTests
    import ..joint_distance_distribution

    import NearestNeighbors
    import StatsBase
    import Distances
    import CausalityToolsBase
    import CausalityToolsBase: causality, CausalityTest, RectangularBinning
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


    # """
    #     causality(source, target, params::CausalityTest)

    # Test whether there is a dynamical influence from `source` to `target` by
    # applying the causality estimator specified by `params`, which must be
    # a valid causality test.

    # ## Examples

    # ```julia
    # x, y = rand(500), rand(500)

    # # Test for a causal relation between `x` and `y` using a cross mapping test
    # # with default parameters
    # causality(x, y, CrossMappingTest())

    # # Test for a causal relation between `x` and `y` using a transfer entropy test
    # # with default parameters over a rectangular binning with five equidistant
    # # intervals along each axis and prediction lag 1
    # binning =  RectangularBinning(5)
    # causality(x, y, TransferOperatorGrid(binning = binning, ηs = 1))

    # # Test for a causal relation between `x` and `y` using a predictive asymmetry
    # # test with a transfer entropy test as the underlying test, using
    # # default parameters over a rectangular binning with ten equidistant
    # # intervals along each axis and prediction lags -10 to 10
    # binning =  RectangularBinning(10)
    # te_test = TransferOperatorGridTest(binning = binning, ηs = -10:10)
    # causality(x, y, PredictiveAsymmetryTest(predictive_test = te_test))
    # ```
    # """
    # function causality end


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

    # Joint distances tests
    # ---------------------
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

    ################################################################
    # Predictive asymmetry causality tests
    ################################################################
    include("tests_predictive_asymmetry/PredictiveAsymmetryTest.jl")

    """
        causality(source, target, test)

    Test for a causal influence from `source` to `target` using the provided causality `test`.

    ## Examples

    ### Real-valued time series

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

    ### Uncertain data

    If your data have uncertainties, use the machinery in `UncertainData.jl` to define uncertain
    values, then use them directly as inputs to `causality`. Uncertain data may be mixed with
    regular real-valued vectors.

    Any combination of real-valued vectors, `Vector{<:AbstractUncertainValue}` and
    `AbstractUncertainValueDataset`s are thus accepted for `source` and `target`. If uncertain data are
    provided, the data are resampled element-wise and the draws are used as input to the causality
    `UncertainIndexValueDataset`s are not yet supported. For more info, see the
    [documentation for `UncertainData.jl`](https://github.com/kahaaga/UncertainData.jl)

    ```julia
    # Vectors of uncertain values. Also, convert to uncertain datasets
    uvals_x = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100]
    uvals_y = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100];
    uvx = UncertainValueDataset(uvals_x)
    uvy = UncertainValueDataset(uvals_y)

    # Draw a single realisation of `uvx` and a single realisation of `uvy`
    x, y = resample(uvx), resample(uvy)

    # A transfer entropy test using the [`TransferOperatorGrid`](@ref) estimator.
    test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1)

    # `causality` accepts all combinations of real-valued and uncertain data types.
    causality(x, y, test_tog)
    causality(x, uvy, test_tog)
    causality(uvx, y, test_tog)
    causality(uvx, uvy, test_tog)
    causality(x, uvals_y, test_tog)
    causality(uvals_x, y, test_tog)
    causality(uvals_x, uvals_y, test_tog)
    ```
    """
    function causality end

    export causality
end



"""
	CausalityTests

A module defining causality tests.
"""
CausalityTests
