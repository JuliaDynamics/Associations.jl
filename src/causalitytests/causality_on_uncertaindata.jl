
import UncertainData: 
    resample,
    AbstractUncertainIndexValueDataset,
    SequentialSamplingConstraint,
    RegularGrid

"""
    causality(source::AbstractUncertainIndexValueDataset, 
        target::AbstractUncertainIndexValueDataset, 
        test::CausalityTest)

Test for a causal influence from `source` to `target` using the provided causality `test`.

The test is performed on a single draw of the values of `source` and `target`, disregarding
the ordering resulting from resampling, but respecting the order of the points in the dataset.

*Note: if the uncertain values furnishing the indices have overlapping supports, you might 
mess up the index-ordering (e.g. time-ordering) of the data points*.

## Example 

Generate some example data series x and y, where x influences y:

```julia 
n_pts = 300
a₁, a₂, b₁, b₂, ξ₁, ξ₂, C₁₂ = 0.7, 0.1, 0.75, 0.2, 0.3, 0.3, 0.5

D = rand(n_pts, 2) 
for t in 5:n_pts
    D[t,1] = a₁*D[t-1,1] - a₂*D[t-4,1] +                ξ₁*rand(Normal(0, 1))
    D[t,2] = b₁*D[t-1,2] - b₂*D[t-4,2] + C₁₂*D[t-1,1] + ξ₂*rand(Normal(0, 1))
end

x, y = D[:, 1], D[:, 2]
```

Add some uncertainties and gather in an `UncertainIndexValueDataset`

```julia
t = 1:n_pts
tu = UncertainValue.(Normal.(t, rand()))
xu = UncertainValue.(Normal.(x, rand()))
yu = UncertainValue.(Normal.(y, rand()))
X = UncertainIndexValueDataset(tu, xu)
Y = UncertainIndexValueDataset(tu, yu)
```

Define a causality test, for example the predictive asymmetry test:

```julia 
# Define causality test 
k, l, m = 1, 1, 1
ηs = -8:8
n_subdivs = floor(Int, n_pts^(1/(k+l+m+1)))
bin = RectangularBinning(n_subdivs)
te_test = VisitationFrequencyTest(k = k, l = l, m = m, binning = bin, ηs = ηs)
pa_test = PredictiveAsymmetryTest(predictive_test = te_test)
```

Run the causality test on a single draw of `X` and a single draw of `Y`:

```
pa_XY = causality(X, Y, pa_test)
pa_YX = causality(Y, X, pa_test)
```

Repeat the test on multiple draws:

```julia
pa_XY = [causality(X, Y, pa_test) for i = 1:100]
pa_YX = [causality(Y, X, pa_test) for i = 1:100]
```
"""
function causality(source::AbstractUncertainIndexValueDataset, 
            target::AbstractUncertainIndexValueDataset, 
            test::CausalityTest)
    
    # @warn """
    # You're running a causality test by resampling `source` and `target`, which both have 
    # uncertain indices, without ensuring proper ordering of the indices. If the 
    # distributions/populations furnishing the indices have overlapping supports, you are 
    # not guaranteed to have the correct ordering, and the results you get might be meaningless! 
    
    # Use the following method instead:
    
    # `causality(source::AbstractUncertainIndexValueDataset, 
    #     target::AbstractUncertainIndexValueDataset, 
    #     constraint::SequentialSamplingConstraint, 
    #     grid::RegularGrid)`. 
    
    # For example, `causality(source, target, StrictlyIncreasing(), RegularGrid(start, stop, step))`
    # will first pick a strictly increasing draw of the age model, interpolate to a regular grid 
    # with spacing `step`, then run the test on the data interpolated to this grid, and finally
    # perform the causality test.
    # """
    idxs_s, vals_s = resample(source)
    idxs_t, vals_t = resample(target)
    causality(vals_s, vals_t, test)
end

"""
    causality(source::AbstractUncertainIndexValueDataset, 
        target::AbstractUncertainIndexValueDataset, 
        test::CausalityTest, constraint::SequentialSamplingConstraint, grid::RegularGrid)

Test for a causal influence from `source` to `target` using the provided causality `test` on 
single draws of `source` and `target` that have been generated according to the provided 
`sequential` sampling constraint. After interpolating both `source` and `target` 
to the provided regular `grid`, the causality test is performed on the interpolated data.

## Example 

First, generate some example data. 

```julia
N = 300
a₁, a₂, b₁, b₂, ξ₁, ξ₂, C₁₂ = 0.7, 0.1, 0.75, 0.2, 0.3, 0.3, 0.5

D = rand(N, 2) 
for t in 5:N
    D[t,1] = a₁*D[t-1,1] - a₂*D[t-3,1] +                ξ₁*rand(Normal(0, 1))
    D[t,2] = b₁*D[t-1,2] - b₂*D[t-2,2] + C₁₂*D[t-1,1] + ξ₂*rand(Normal(0, 1))
end
```

Gather time series and add some uncertainties to them:

```julia
ts = collect(1:N)
x, y = D[:, 1], D[:, 2]
t = UncertainValue.(Normal.(ts .+ 0.2 .* rand(N), rand(N)))
uvalx = UncertainValue.(Normal.(x, rand()))
uvaly = UncertainValue.(Normal.(y, rand()))

X = UncertainIndexValueDataset(t, uvalx)
Y = UncertainIndexValueDataset(t, uvaly)
```

Now, resample both the indices and values of `X` and `Y` in a strictly 
increasing manner according to the indices, interpolate both 
to the same regular grid, and perform a causality test.

```julia
n_bins = ceil(Int, N^(1/4))
te_test = VisitationFrequencyTest(binning = RectangularBinning(n_bins), ηs = -5:5)
pa_test = PredictiveAsymmetryTest(predictive_test = te_test)

# Sample a strictly increasing realisation of the indices, interpolate,
# then perform the causality test on the interpolated data.
causality(X, Y, pa_test, StrictlyIncreasing(), RegularGrid(1:1:N))
```
"""
function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest,
    constraint::SequentialSamplingConstraint,
    grid::RegularGrid)
    
    idxs_s, vals_s = resample(source, constraint, grid)
    idxs_t, vals_t = resample(target, constraint, grid)
    causality(vals_s, vals_t, test)
end

"""
    causality(x, y, test::CausalityTest)

Test for a causal influence from `source` to `target` using the provided causality `test`.

Both `x` and `y` can be a real-valued vector, `Vector{<:AbstractUncertainValue}` 
or a `AbstractUncertainValueDataset`. If either `x`, `y` or both are uncertain, then the test 
is applied to a single draw from the uncertain data. 

## Examples 

Generate some example data series x and y, where x influences y:

```julia
n_pts = 300
a₁, a₂, b₁, b₂, ξ₁, ξ₂, C₁₂ = 0.7, 0.1, 0.75, 0.2, 0.3, 0.3, 0.5

D = rand(n_pts, 2) 
for t in 5:n_pts
    D[t,1] = a₁*D[t-1,1] - a₂*D[t-3,1] +                ξ₁*rand(Normal(0, 1))
    D[t,2] = b₁*D[t-1,2] - b₂*D[t-2,2] + C₁₂*D[t-1,1] + ξ₂*rand(Normal(0, 1))
end
```

Gather time series and add some uncertainties to them:

```julia
x, y = D[:, 1], D[:, 2]
uvalx = UncertainValue.(Normal.(x, rand()))
uvaly = UncertainValue.(Normal.(y, rand()))

xd = UncertainValueDataset(uvalx)
yd = UncertainValueDataset(uvaly)
```

Any combination of certain and uncertain values will work:

```julia
causality(x, y, pa_test)
causality(x, yd, pa_test)
causality(xd, yd, pa_test)
causality(x, uvaly, pa_test)
causality(uvalx, uvaly, pa_test)
```

On multiple realisations of the uncertain `yd`, but fixing `x`:

```julia
[causality(x, yd, pa_test) for i = 1:100]
```
"""
function causality(source::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset},
        target::Union{AbstractVector, Vector{<:AbstractUncertainValue}, AbstractUncertainValueDataset},
        test::CausalityTest)
end
 
function causality(source::AbstractVector{T}, 
    target::Vector{<:AbstractUncertainValue}, 
    test::CausalityTest) where T

    causality(source, resample(target), test)
end

function causality(source::Vector{<:AbstractUncertainValue}, 
    target::AbstractVector{T}, 
    test::CausalityTest) where T

    causality(resample(source), target, test)
end

function causality(source::AbstractVector{T}, 
    target::AbstractUncertainValueDataset, 
    test::CausalityTest) where T

    causality(source, resample(target), test)
end

function causality(source::AbstractUncertainValueDataset, 
    target::AbstractVector{T}, 
    test::CausalityTest) where T

    causality(resample(source), target, test)
end