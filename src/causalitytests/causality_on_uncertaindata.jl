
import UncertainData: 
    resample,
    AbstractUncertainIndexValueDataset,
    SequentialSamplingConstraint,
    RegularGrid

function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest)
    
    @warn """
    You're running a causality test by resampling `source` and `target`, which both have 
    uncertain indices, without ensuring proper ordering of the indices. If the 
    distributions/populations furnishing the indices have overlapping supports, you are 
    not guaranteed to have the correct ordering, and the results you get might be meaningless! 
    
    Use the following method instead:
    
    `causality(source::AbstractUncertainIndexValueDataset, 
        target::AbstractUncertainIndexValueDataset, 
        constraint::SequentialSamplingConstraint, 
        grid::RegularGrid)`. 
    
    For example, `causality(source, target, StrictlyIncreasing(), RegularGrid(start, stop, step))`
    will first pick a strictly increasing draw of the age model, interpolate to a regular grid 
    with spacing `step`, then run the test on the data interpolated to this grid, and finally
    perform the causality test.
    """
    idxs_s, vals_s = resample(source)
    idxs_t, vals_t = resample(target)
    causality(vals_s, vals_t, test)
end


"""
    causality(source::AbstractUncertainIndexValueDataset, 
        target::AbstractUncertainIndexValueDataset, 
        test::CausalityTest,
        constraint::SequentialSamplingConstraint,
        grid::RegularGrid)

Test for a causal influence from `source` to `target` using the provided causality `test`.
Runs the test on a single realisation pair of `source` and `target` using the provided 
sequential sampling `constraint`, interpolates both `source` and `target` to the provided 
regular `grid`, then performs the causality test on the interpolated data.

## Examples

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
function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest,
    constraint::SCT,
    grid::RegularGrid) where SCT <: SequentialSamplingConstraint
    
    idxs_s, vals_s = resample(source, constraint, grid)
    idxs_t, vals_t = resample(target, constraint, grid)
    causality(vals_s, vals_t, test)
end