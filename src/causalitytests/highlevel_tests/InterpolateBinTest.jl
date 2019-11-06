import UncertainData: InterpolateAndBin

"""
    InterpolateBinTest(test::CausalityTest, interpolate_bin_resampling::InterpolateAndBin)
    InterpolateBinTest(test::CausalityTest, interpolate_bin_resampling::InterpolateAndBin, n::Int)

A causality test where the data is interpolated and then binned before the test is applied.

`n` controls the number of independent resamplings over which the test 
is performed (if not provided, `n` is set to `1` by default).

## Examples 

Record `N` points from the built-in `ar1_unidir` system, collect the 1st and 2nd
variables as `X` and `Y` and add some uncertainties to both the indices and the values.

```julia 
N = 300
sys = ar1_unidir(c_xy = 0.8)
X, Y = example_uncertain_indexvalue_datasets(sys, N, (1, 2),
    d_xval = Uniform(0.001, 0.05), d_yval = Uniform(0.001, 0.05));
```

Say we want to interpolate both `X` and `Y` to some very fine grid 
using linear interpolation, then bin to a coarser grid and summarise
each bin using the mean of the values in each bin.

```julia
# Define interpolation grid over the range of available index values
tmin = max(minimum(mean.(X.indices)), minimum(mean.(Y.indices)))
tmax = max(maximum(mean.(X.indices)), maximum(mean.(Y.indices)))
intp_grid = tmin:0.01:tmax

# Define binning grid
left_bin_edges = tmin:5:tmax

# Define the InterpolateAndBin instance
intp_bin = InterpolateAndBin(mean, left_bin_edges, Linear(), intp_grid, Flat(OnGrid()))
```
"""
struct InterpolateBinTest{
        CT <: CausalityTest, 
        IB <: InterpolateAndBin{L} where L} <: PreprocessedDataCausalityTest
    test::CT
    interpolate_bin_resampling::IB
    n::Int
end

function InterpolateBinTest(test, interpolate_bin_resampling)
    InterpolateBinTest(test, interpolate_bin_resampling, 1)
end;

export InterpolateBinTest