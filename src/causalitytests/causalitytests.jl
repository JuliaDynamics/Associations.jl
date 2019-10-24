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


    """
        causality(source, target, test::CausalityTest)
        causality(source::AbstractUncertainIndexValueDataset, 
            target::AbstractUncertainIndexValueDataset, 
            test::CausalityTest, 
            constraint::SequentialSamplingConstraint, 
            grid::RegularGrid)

    Test for a causal influence from `source` to `target` using the provided causality `test`.
    
    # Uncertainty handling 

    All high-level causality test are integrated with the uncertainty handling machinery 
    in [UncertainData.jl](https://github.com/kahaaga/UncertainData.jl). See the
    [UncertainData.jl documentation](https://kahaaga.github.io/UncertainData.jl/dev/)
    for more details.

    # Valid input types for `source` and `target`
    
    # Uncertainties in values 
    
    `source` and `target` may be any combination of real-valued vectors, 
    `Vector{<:AbstractUncertainValue}`, or `AbstractUncertainValueDataset`. If either 
    input is uncertain, then the causality test is performed on a random draw of 
    that uncertain dataset. 

    ## Uncertainties in indices and values

    If `source` and `target` are *both* instances of some subtype of 
    `UncertainData.AbstractUncertainIndexValueDataset`, then the test is performed 
    on a single realisation pair of `source` and `target`, drawn using the provided 
    sequential sampling `constraint`. After interpolating both `source` and `target` 
    to the provided regular `grid`, the causality test is done on the interpolated data.
    
    # Examples

    ## Real-valued time series

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

    # Uncertain data
    
    ## Only values are uncertain 

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

    ## Indices and values are uncertain 
    
    ### A simple example

    ```julia
    using CausalityTools, UncertainData

    # Define a system and generate some time series
    system = logistic2_unidir(c_xy = 0.5, r₁ = 3.7, r₂ = 3.8)
    orbit = trajectory(system, 150, Ttr = 1000)

    # Add uncertainties to the time series
    x_uncertain = [UncertainValue(Normal, x, rand(Uniform(0.001, 0.07))) for x in orbit[:, 1]]
    y_uncertain = [UncertainValue(Normal, y, rand(Uniform(0.001, 0.07))) for y in orbit[:, 2]]
    x = UncertainValueDataset(x_uncertain)
    y = UncertainValueDataset(y_uncertain)

    # Add uncertainties to the time indices 
    time = [UncertainValue(Normal, i, rand()) for i = 1:length(x)];
    timeinds = UncertainIndexDataset(time)

    # Define some `UncertainIndexValueDataset`s
    X = UncertainIndexValueDataset(timeinds, x)
    Y = UncertainIndexValueDataset(timeinds, y);

    # Plot the data
    pX = plot(lw = 0.5, c = :blue, size = (1000, 300), xlabel = "Value", ylabel = "X", label = "")
    plot!(pX, mean.(X.indices), mean.(X.values), lw = 2, c = :blue)
    plot!(pX, X, c = :blue)

    pY = plot(lw = 0.5, c = :red, size = (1000, 300), xlabel = "Value", ylabel = "Y", label = "")
    plot!(pY, mean.(Y.indices), mean.(Y.values), lw = 2, c = :red)
    plot!(pY, Y, c = :red)

    plot(pX, pY, layout = (2, 1), legend = false)

    # We have 101 values, so let's interpolate the random draw to 
    # back to the grid given by the sequence 0:1:101
    grid = RegularGrid(0, 101, 1)

    # Define causality test 
    ηs = -5:5
    te_test = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = ηs)
    test = PredictiveAsymmetryTest(predictive_test = te_test)

    # Sample a strictly increasing realisation of the indices, interpolate,
    # then perform the causality test on the interpolated data.
    causality(X, Y, test, StrictlyIncreasing(), RegularGrid(0, 101, 1))
    ```

    This single draw might result in either positive or negative values, but that's 
    not really informative. When we have access to the uncertainties in the data 
    indices and/or values, we should use them to obtain uncertainties on our 
    causality statistic. 

    In the following example, we resample from the age-value models for both datasets 
    and compute the predictive asymmetry over an ensemble of realisations that are 
    drawn from within the precision of the data, imposing strictly increasing values 
    for the age model.

    ### A more complex example: multiple realisations

    ```julia
    using CausalityTools, UncertainData, Plots 

    # Define a system and generate some time series
    system = logistic2_unidir(c_xy = 0.7, r₁ = 3.7, r₂ = 3.8)
    n_points =  500
    orbit = trajectory(system, n_points, Ttr = 1000)

    # Add uncertainties to the time series
    x_uncertain = [UncertainValue(Normal, x, rand(Uniform(0.001, 0.07))) for x in orbit[:, 1]]
    y_uncertain = [UncertainValue(Normal, y, rand(Uniform(0.001, 0.07))) for y in orbit[:, 2]]
    x = UncertainValueDataset(x_uncertain)
    y = UncertainValueDataset(y_uncertain)

    # Add uncertainties to the time indices 
    time = [UncertainValue(Normal, i, rand()) for i = 1:length(x)];
    timeinds = UncertainIndexDataset(time)

    # Define some `UncertainIndexValueDataset`s
    X = UncertainIndexValueDataset(timeinds, x)
    Y = UncertainIndexValueDataset(timeinds, y);

    # We have 101 values, so let's interpolate the random draw to 
    # back to the grid given by the sequence 0:1:101
    grid = RegularGrid(0, 101, 1)

    # Define causality test 
    η_max = 5
    ηs = -η_max:η_max
    n_bins = ceil(Int, n_points^(1/4))

    te_test = VisitationFrequencyTest(binning = RectangularBinning(n_bins), ηs = ηs)

    test = PredictiveAsymmetryTest(predictive_test = te_test)

    # Sample a strictly increasing realisation of the indices, interpolate,
    # then perform the causality test on the interpolated data.
    causality(X, Y, test, StrictlyIncreasing(), RegularGrid(0, 101, 1))

    # Perform the test on 300 realisations of the age-data uncertainty model 
    tstep = 1
    grid = RegularGrid(0, length(X), tstep)
    causality_xtoy = [causality(X, Y, test, StrictlyIncreasing(), grid) for i = 1:100]
    causality_ytox = [causality(Y, X, test, StrictlyIncreasing(), grid) for i = 1:100]

    # Collect the result vectors and summarise (mean and standard deviation) for 
    # each prediction lag
    xtoys = hcat(causality_xtoy...,)
    ytoxs = hcat(causality_ytox...,);

    xtoys_mean = dropdims(mean(xtoys, dims = 2), dims = 2)
    ytoxs_mean = dropdims(mean(ytoxs, dims = 2), dims = 2)
    xtoys_std = dropdims(std(xtoys, dims = 2), dims = 2)
    ytoxs_std = dropdims(std(ytoxs, dims = 2), dims = 2);


    # Plot the results, so that we can interpret them. Remember that the 
    # prediction lags are in units of the interpolation grid.
    p = plot(xlabel = "Prediction lag", ylabel = "Predictive asymmetry (bits)",
    xlims = (1, η_max))
    plot!(p, 1:η_max, xtoys_mean, ribbons = xtoys_std, c = :black, label = "x to y")
    plot!(p, 1:η_max, ytoxs_mean, ribbons = ytoxs_std, c = :red, label = "y to x")
    plot!(p, 1:η_max, xtoys_mean, lw = 2, c = :black, label = "")
    plot!(p, 1:η_max, ytoxs_mean, lw = 2, c = :red, label = "")
    hline!([0], ls = :dash, lw = 2, lc = :grey, label = "") # zero line
    ```

    ### When only one dataset has uncertain time indices 

    When only one time series has uncertainties in the time indices, we still
    construct a `UncertainValueIndexDataset`, but represent the indices 
    as `CertainValue` instances.

    ```julia
    using CausalityTools, UncertainData, Plots 

    # Define a system and generate some time series
    system = logistic2_unidir(c_xy = 0.7, r₁ = 3.66, r₂ = 3.81)
    n_points =  400
    orbit = trajectory(system, n_points, Ttr = 1000)

    # Add uncertainties to the time series values
    x_uncertain = [UncertainValue(Normal, x, rand(Uniform(0.001, 0.07))) for x in orbit[:, 1]]
    y_uncertain = [UncertainValue(Normal, y, rand(Uniform(0.001, 0.07))) for y in orbit[:, 2]]
    x = UncertainValueDataset(x_uncertain)
    y = UncertainValueDataset(y_uncertain)

    # Add uncertainties to only the time indices for. Take the 
    # time indices for y as certain.
    time_uncertain = [UncertainValue(Normal, i, 0.5) for i = 1:length(x)];
    time_certain = [CertainValue(i) for i = 1:length(x)];
    timeinds_x = UncertainIndexDataset(time_uncertain)
    timeinds_y = UncertainIndexDataset(time_certain)

    X = UncertainIndexValueDataset(timeinds_x, x)
    Y = UncertainIndexValueDataset(timeinds_y, y)

    # We have 101 values, so let's interpolate the random draw to 
    # back to the grid given by the sequence 0:1:101
    grid = RegularGrid(0, 101, 1)

    # Define causality test 
    η_max = 5
    ηs = -η_max:η_max
    n_bins = ceil(Int, n_points^(1/4))

    te_test = TransferOperatorGridTest(binning = RectangularBinning(n_bins), ηs = ηs)

    test = PredictiveAsymmetryTest(predictive_test = te_test)

    # Sample strictly increasing realisations of the indices, interpolate,
    # then perform the causality test on the interpolated data.
    tstep = 1
    grid = RegularGrid(0, length(X), tstep)
    causality_xtoy = [causality(X, Y, test, StrictlyIncreasing(), grid) for i = 1:50]
    causality_ytox = [causality(Y, X, test, StrictlyIncreasing(), grid) for i = 1:50]

    # Collect the result vectors and summarise (mean and standard deviation) for 
    # each prediction lag
    xtoys = hcat(causality_xtoy...,)
    ytoxs = hcat(causality_ytox...,);

    xtoys_mean = dropdims(mean(xtoys, dims = 2), dims = 2)
    ytoxs_mean = dropdims(mean(ytoxs, dims = 2), dims = 2)
    xtoys_std = dropdims(std(xtoys, dims = 2), dims = 2)
    ytoxs_std = dropdims(std(ytoxs, dims = 2), dims = 2);

    # Plot the results, so that we can interpret them
    p = plot(xlabel = "Prediction lag", ylabel = "Predictive asymmetry (bits)",
    xlims = (1, η_max))
    plot!(p, 1:η_max, xtoys_mean, ribbons = xtoys_std, c = :black, label = "x to y")
    plot!(p, 1:η_max, ytoxs_mean, ribbons = ytoxs_std, c = :red, label = "y to x")
    plot!(p, 1:η_max, xtoys_mean, lw = 2, c = :black, label = "")
    plot!(p, 1:η_max, ytoxs_mean, lw = 2, c = :red, label = "")
    hline!([0], ls = :dash, lw = 2, lc = :grey, label = "") # zero line
    ```
    """
    function causality end


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
    include("PreprocessedDataCausalityTest.jl")
    include("BinnedDataCausalityTest.jl")

    export causality
end



"""
	CausalityTests

A module defining causality tests.
"""
CausalityTests
