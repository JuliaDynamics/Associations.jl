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

