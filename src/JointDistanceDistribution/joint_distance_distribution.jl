using HypothesisTests
using Distances
using DelayEmbeddings

function normalise_minmax(x, vmin, vmax)
    (x - vmin)/(vmax - vmin)
end

"""
    jdd(source, target; distance_metric = SqEuclidean(), 
        B::Int = 10, D::Int = 2, τ::Int = 1) → Vector{Float64}

Compute the joint distance distribution [1] from `source` to `target` using 
the provided `distance_metric`, with `B` controlling the number of subintervals, 
`D` the embedding dimension and `τ` the embedding lag.

## Example

```julia
using CausalityTools
x, y = rand(1000), rand(1000)

jdd(x, y)
```

## Keyword arguments

- **`distance_metric::Metric`**: An instance of a valid distance metric from `Distances.jl`. 
    Defaults to `SqEuclidean()`.
- **`B`**: The number of equidistant subintervals to divide the interval `[0, 1]` into
    when comparing the normalised distances. 
- **`D`**: Embedding dimension.
- **`τ`**: Embedding delay.


## References 
[1] Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of Nonlinear Science 28.7 (2018): 075302.
"""
function jdd(source, target;
        distance_metric = Euclidean(),
        B::Int = 10, 
        D::Int = 2, 
        τ::Int = 1)
    
    js = ([1 for i = 1:D]...,)
    τs = (collect(0:-τ:-(D-1)*τ)...,)
    Ex = DelayEmbeddings.genembed(source, τs, js)
    Ey = DelayEmbeddings.genembed(target, τs, js)
    Mx = DelayEmbeddings.Matrix(Ex)
    My = DelayEmbeddings.Matrix(Ey)
    
    npts = length(Ex)
    Dx = pairwise(distance_metric, Mx, Mx, dims = 1)
    Dy = pairwise(distance_metric, My, My, dims = 1)
    
    # Normalise the distances to the interval [0, 1]
    Dx_min = minimum(Dx[Dx .> 0])
    Dy_min = minimum(Dy[Dy .> 0])
    Dx_max = maximum(Dx[Dx .> 0])
    Dy_max = maximum(Dy[Dy .> 0])
    
    Dx_norm = zeros(Float64, size(Dx))
    Dy_norm = zeros(Float64, size(Dy))

    for i in LinearIndices(Dx[Dx .> 0])
        Dx_norm[i] = normalise_minmax(Dx[i], Dx_min, Dx_max)
    end

    for i in LinearIndices(Dy[Dy .> 0])
        Dy_norm[i] = normalise_minmax(Dy[i], Dy_min, Dy_max)
    end
    
    mins_δ_yi_yj = fill(2.0, 2*B)

    for (k, b) in enumerate(1:2*B)
        bmin = (b-1)/(2*B)
        bmax = b/(2*B)
        
        # Find the indices of all pairs (i, j) in Dx whose distances fall inside the interval Ib
        #idxs_Dxs_in_Ib = findall(Dx[])

        # We don't need to store any indices or distances explicitly, but only 
        # keep track of whether a smaller distance has has been detected. 
        # The maximum possible distance after normalisation is 1.0, so this 
        # value can only decrease as we update. 
        min_δ_yi_yj = 1.0

        for i = 1:npts
            for j = (i+1):npts
                δ_xi_xj = Dx_norm[i, j]
                if bmin < δ_xi_xj <= bmax
                    δ_yi_yj = Dy_norm[i, j]
                    if δ_yi_yj < min_δ_yi_yj
                        min_δ_yi_yj = δ_yi_yj
                    end
                end
            end
        end
        mins_δ_yi_yj[k] = min_δ_yi_yj
    end
    
    Δjdd = [mins_δ_yi_yj[B + i] - mins_δ_yi_yj[i] for i in 1:B]
    
    return Δjdd
end

"""
    jdd(test::OneSampleTTest, source, target;
        distance_metric = SqEuclidean(), B::Int = 10, D::Int = 2, τ::Int = 1, 
        μ0 = 0.0) → OneSampleTTest

Perform a one sample t-test to check that the joint distance distribution [1] 
computed from `source` to `target` is biased towards positive values, using the null 
hypothesis that the mean of the distribution is `μ0`.

The interpretation of the t-test is that if we can reject the null, then the 
joint distance distribution is biased towards positive values, and then there 
exists an underlying coupling from `source` to `target`. 


## Example

```julia 
using CausalityTools, HypothesisTests
x, y = rand(1000), rand(1000)

jdd(OneSampleTTest, x, y)
```

which gives 

```julia
One sample t-test
-----------------
Population details:
    parameter of interest:   Mean
    value under h_0:         0.0
    point estimate:          0.06361857324022721
    95% confidence interval: (0.0185, 0.1087)

Test summary:
    outcome with 95% confidence: reject h_0
    two-sided p-value:           0.0082

Details:
    number of observations:   20
    t-statistic:              2.9517208721082873
    degrees of freedom:       19
    empirical standard error: 0.0215530451545668
```


The lower bound of the confidence interval for the mean of the joint 
distance distribution is `0.0185` at confidence level `α = 0.05`. The 
meaning that the test falsely detected causality from `x` to `y`
between these two random time series. To get the confidence intervals
at confidence level `α`, use `confinf(jdd, α)`. If you just want the 
p-value at 95% confidence, use `pvalue(jdd, tail = :left)`

## Keyword arguments

- **`distance_metric::Metric`**: An instance of a valid distance metric from `Distances.jl`. 
    Defaults to `SqEuclidean()`.
- **`B`**: The number of equidistant subintervals to divide the interval `[0, 1]` into
    when comparing the normalised distances. 
- **`D`**: Embedding dimension.
- **`τ`**: Embedding delay.
- **`μ0`**: The hypothetical mean value of the joint distance distribution if there 
    is no coupling between `x` and `y` (default is `μ0 = 0.0`).

## References 
[1] Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of Nonlinear Science 28.7 (2018): 075302.
"""
function jdd(test::Type{OneSampleTTest}, source, target; 
        distance_metric = SqEuclidean(), 
        B::Int = 10, 
        D::Int = 2, τ::Int = 1, 
        μ0 = 0.0)

    Δjdd = jdd(source, target, 
        distance_metric = distance_metric, 
        B = B, 
        D = D, τ = τ)
    
    OneSampleTTest(Δjdd, μ0)
end 

export jdd