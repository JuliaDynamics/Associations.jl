import UncertainData: 
    AbstractBinnedResampling, 
    AbstractBinnedUncertainValueResampling,
    AbstractBinnedSummarisedResampling,
    AbstractUncertainIndexValueDataset,
    resample

export BinnedDataCausalityTest, causality

"""
    BinnedDataCausalityTest(test::CT, binning::AbstractBinnedUncertainValueResampling, n_realizations::Int)
    BinnedDataCausalityTest(test::CT, binning::AbstractBinnedSummarisedResampling)

A causality test where the data are binned before applying the test. If `binning` is some 
binning scheme that returns an uncertain value for each 
bin (e.g. `BinnedResampling` or `BinnedWeightedResampling`), then the `test` is 
applied `n_realizations` times. If `binning` returns a single value for each bin 
(e.g. `BinnedMeanResampling` or `BinnedMeanWeightedResampling`, the `test` is applied
only once.

## Fields

- **`test::CausalityTest`**. An instance of a causality test, e.g. `VisitationFrequencyTest`, 
    `PredictiveAsymmetryTest`, `CrossMappingTest` or `JointDistanceDistributionTest`.

- **`binning::AbstractBinnedResampling`**. An instance of a resampling scheme 
    indicating the type of binning, e.g `BinnedResampling`, `BinnedWeightedResampling`, 
    `BinnedMeanResampling` or `BinnedMeanWeightedResampling`.

- **`n_realizations::Int`**: The number of independent draws of the binned data set over
    which to compute the causality `test`. Only taken into consideration of `binning` is 
    a scheme resulting in each bin being represented by an uncertain value. If bins are 
    summarized, then the `test` is always applied only once.

## Examples 

### Binned convergent cross mapping test

Let's say we want to bin the data by drawing n_draws = 5000 realiastions of 
each uncertain data point, then assign the draws to the correct bins, and
finally get a kernel density estimate to the distribution of values in each bin.

```julia 
grid = 0:10:1000 # left bin edges
n_draw = 5000 # sample each point 5000 times and distribute among bins
binning = BinnedResampling(grid, 5000)

ccm_test = ConvergentCrossMappingTest(timeseries_lengths = 15:5:100 |> collect, 
    n_reps = 100, libsize = 100, τ = 1, replace = true)

test = BinnedDataCausalityTest(ccm_test, binning, n_realizations)
```

### Binned transfer entropy test

Let's say we want to bin the data by drawing n_draws = 7500 realiastions of 
each uncertain data point, then assign the draws to the correct bins, and
finally get a kernel density estimate to the distribution of values in each bin.

```julia 
grid = 0:5:1000 # left bin edges
n_draws = 7500
binning = BinnedResampling(grid, n_draws)

# A transfer entropy causality test using the visitation frequency estimator
state_space_binning = RectangularBinning(5)
te_test = VisitationFrequencyTest(binning = state_space_binning, ηs = 1)

test = BinnedDataCausalityTest(te_test, binning, n_realizations)
```
"""
struct BinnedDataCausalityTest{
    CT <: CausalityTest, 
    BR <: AbstractBinnedResampling} <: PreprocessedDataCausalityTest

    test::CT
    binning::BR
    n_realizations::Int
end

function BinnedDataCausalityTest(test::CR, binning::BR) where {
    CR <: CausalityTest, BR <: AbstractBinnedSummarisedResampling}

    BinnedDataCausalityTest(test, binning, 1) 
end



"""
    causality(x::AbstractUncertainIndexValueDataset,
        y::AbstractUncertainIndexValueDataset,
        test::BinnedDataCausalityTest)

Apply a causality test to `x` and `y`, which are both data series 
with uncertainties in both indices and values. To get the data 
on an equally-spaced temporal grid, the data are first binned
according to the instructions in `test.binning`. 

## Binning methods 

- If `test.binning` results in an uncertain value for each bin, then the causality 
    test is applied `test.n_realizations` times. Examples are [`BinnedResampling`](@ref) 
    and [`BinnedWeightedResampling`](@ref), which both return a KDE estimate to the 
    distribution of values in each bin.

- If `test.binning` returns a summary statistic for each bin, then the causality test is 
    applied to the summarised time series. Examples are [`BinnedMeanResampling`](@ref) 
    and [`BinnedMeanWeightedResampling`](@ref), which both return the mean of each bin.

## Example

### Some example data 

Let's say we have the following datasets `X` and `Y`: 

```julia
using CausalityTools, UncertainData 

sys = ar1_unidir(uᵢ = [0.1, 0.1], c_xy = 0.41)
vars = (1, 2) # ar1_unidir has only two variables, X and Y
n_steps = 100 # the number of points in the time series
tstep = 10 # the mean of each time value is stepped by `tstep`

X, Y = example_uncertain_indexvalue_datasets(sys, n_steps, vars, tstep = 10, 
    d_xind = Uniform(7.5, 15.5), 
    d_yind = Uniform(5.5, 15.5), 
    d_xval = Uniform(0.1, 0.5));
```

Now we can perform the causality tests on our uncertain data, either over 
the bin means, or over an ensemble of independent realisations of 
the dataset given the distribution of points in each bin.

### Summary statistic for each bin

#### PredictiveAsymmetryTest

If the binning method returns a summary statistic (e.g. the mean) for each bin, 
then the causality test is applied exactly once to the bin summaries.

```julia
# Bin the data by drawing n_draws = 5000 realizations of each uncertain data 
# point, then assign the draws to the correct bins and get a distribution 
# of values in each bin.
binning = BinnedResampling(0:10:1000, 5000)

# A predictive asymmetry test. For the predictions, we'll use a transfer entropy 
# test, which uses the visitation frequency estimator. We'll predict 5 time steps 
# forwards and backwards in time (ηs = -5:5). For the state space binning, we'll
# determine how many intervals to split each axis into from the number of points
# in the time series.
k, l, m = 1, 1, 1 # embedding parameters
n_subdivisions = floor(Int, length(0:10:1000)^(1/(k + l + m + 1)))
state_space_binning = RectangularBinning(n_subdivisions)
ηs = -5:5
te_test = VisitationFrequencyTest(k = k, l = l, m = m,
    binning = state_space_binning, 
    ηs = ηs, b = 2) # use base-2 logarithms
pa_test = PredictiveAsymmetryTest(predictive_test = te_test)

# Compute transfer entropy in both directions for the bin means.
# We still have to specify the number of realisations for the test, 
# but it is ignored, so it doesn't matter what number we input here.
test = BinnedDataCausalityTest(pa_test, binning, 1)

# Perform the tests
tes_xy = causality(X, Y, test)
tes_yx = causality(Y, X, test)

# Plot the results
plot(xlabel = "Prediction lag (eta)", ylabel = "Predictive asymmetry (bits)")
plot!(ηs[ηs .> 0], tes_xy, label = "x -> y (binned)", c = :black)
plot!(ηs[ηs .> 0], tes_yx, label = "y -> x (binned)", c = :red)
hline!([0], lw = 2, ls = :dot, α = 0.5, label = "", c = :grey)
```

### Uncertain values representing each bin


If the binning method returns an uncertain value for each bin, then the 
causality test is applied to `test.n_realizations` independent draws 
of the binned dataset.

#### Predictive asymmetry test

```julia
# Bin the data by drawing n_draws = 5000 realizations of each uncertain data 
# point, then assign the draws to the correct bins and get a distribution 
# of values in each bin.
grid = 0:5:1000
binning = BinnedResampling(grid, 5000)

# A predictive asymmetry test. For the predictions, we'll use a transfer entropy 
# test, which uses the visitation frequency estimator. We'll predict 5 time steps 
# forwards and backwards in time (ηs = -5:5). For the state space binning, we'll
# determine how many intervals to split each axis into from the number of points
# in the time series.
k, l, m = 1, 1, 1 # embedding parameters
n_subdivisions = floor(Int, length(0:10:1000)^(1/(k + l + m + 1)))
state_space_binning = RectangularBinning(n_subdivisions)
ηs = -5:5
te_test = VisitationFrequencyTest(k = k, l = l, m = m,
    binning = state_space_binning, 
    ηs = ηs, b = 2) # use base-2 logarithms
pa_test = PredictiveAsymmetryTest(predictive_test = te_test)

# Compute transfer entropy in both directions over 50 independent realizations 
# of the binned dataset.
n_realizations = 50
test = BinnedDataCausalityTest(pa_test, binning, n_realizations)

# Perform the tests on the bin means. `tes_xy` now contains 50 independent 
# computation of the predictive asymmetry at predictive lags 1 to 5, 
# and `tes_yx` contains the same but for the other direction
tes_xy = causality(X, Y, test)
tes_yx = causality(Y, X, test)

# Gather results in a matrix and compute means and standard deviations 
# for the predictive asymmetries at each prediction lag
M_xy = hcat(tes_xy...,)
M_yx = hcat(tes_yx...,)

means_xy = mean(M_xy, dims = 2)[:, 1]
means_yx = mean(M_yx, dims = 2)[:, 1]
stdevs_xy = std(M_yx, dims = 2)[:, 1]
stdevs_yx = std(M_yx, dims = 2)[:, 1]

# Plot the predictive asymmetry as a function of prediction lag
plot(xlabel = "Prediction lag (eta)", ylabel = "Predictive asymmetry (bits)")
plot!(ηs[ηs .> 0], means_xy, ribbon = stdevs_xy, label = "x -> y (binned)", c = :black)
plot!(ηs[ηs .> 0], means_yx, ribbon = stdevs_yx, label = "y -> x (binned)", c = :red)
hline!([0], lw = 2, ls = :dot, α = 0.5, label = "", c = :grey)
```

#### Convergent cross mapping test

```julia
using Plots 

# Bin the data by drawing n_draws = 5000 realizations of each uncertain data 
# point, then assign the draws to the correct bins and get a distribution 
# of values in each bin.
binning = BinnedResampling(0:10:1000, 5000)

# A convergent cross mapping causality test
ts_lengths = 15:5:100 |> collect
causalitytest = ConvergentCrossMappingTest(
    timeseries_lengths = ts_lengths, 
    n_reps = 100, 
    libsize = 100, 
    τ = 1, 
    replace = true)

# Run ccm in both directions over 50 independent realizations of the binned datasets
n_realizations = 50
test = BinnedDataCausalityTest(causalitytest, binning, n_realizations)

ccms_xy = causality(X, Y, test);
ccms_yx = causality(Y, X, test);

realization_means_xy = [mean.(realization) for realization in ccms_xy]
realization_means_yx = [mean.(realization) for realization in ccms_yx]
M_xy = hcat(realization_means_xy...,)
M_yx = hcat(realization_means_yx...,)

plot(xlabel = "Time series length", ylabel = "Cross map skill")
plot!(ts_lengths, mean(M_xy, dims = 2)[:, 1], 
    ribbon = std(M_xy, dims = 2)[:, 1], 
    label = "x -> y (binned)")
plot!(ts_lengths, mean(M_yx, dims = 2)[:, 1], 
    ribbon = std(M_yx, dims = 2)[:, 1], 
    label = "y -> x (binned)")
```
"""
function causality(x::AbstractUncertainIndexValueDataset,
    y::AbstractUncertainIndexValueDataset,
    test::BinnedDataCausalityTest)
end


function causality(x::AbstractUncertainIndexValueDataset,
    y::AbstractUncertainIndexValueDataset,
    test::BinnedDataCausalityTest{CT, BR}) where {
        CT <: CausalityTest, BR <: AbstractBinnedUncertainValueResampling}

    binned_x = resample(x, test.binning)
    binned_y = resample(y, test.binning)
    
    [causality(resample(binned_x.values), resample(binned_y.values), 
            test.test) for i = 1:test.n_realizations]
end

function causality(x::AbstractUncertainIndexValueDataset,
    y::AbstractUncertainIndexValueDataset,
    test::BinnedDataCausalityTest{CT, BR}) where {
        CT <: CausalityTest, BR <: AbstractBinnedSummarisedResampling}

    binned_and_summarised_x = resample(x, test.binning)
    binned_and_summarised_y = resample(y, test.binning)
    
    causality(binned_and_summarised_x, binned_and_summarised_y, test.test) 
end
