import Distances: Metric, SqEuclidean, Euclidean

"""
    SMeasureTest(m, τ, K)
    SMeasureTest(m::Int, τ, K; metric::Metric = SqEuclidean(), 
        tree_metric::Metric = Euclidean())

S-measure test [1] for the directional dependence between data series.

## Algorithm

1. Create an `m-`dimensional embeddings of both `x` and `y`, resulting in 
    `N` different `m`-dimensional embedding points
    ``X = \\{x_1, x_2, \\ldots, x_N \\}`` and ``X = \\{y_1, y_2, \\ldots, y_N \\}``.
    `τ` controls the embedding lag.
2. Let ``r_{i,j}`` and ``s_{i,j}`` be the indices of the `k`-th nearest neighbors 
    of ``x_i`` and ``y_i``, respectively.
3. Compute the the mean squared Euclidean distance to the ``k`` nearest neighbors 
    for each ``x_i``, using the indices ``r_{i, j}``.

```math
R_i^{(k)}(x) = \\dfrac{1}{k} \\sum_{i=1}^{k}(x_i, x_{r_{i,j}})^2
```

- Compute the y-conditioned mean squared Euclidean distance to the ``k`` nearest 
    neighbors for each ``x_i``, now using the indices ``s_{i,j}``.

```math
R_i^{(k)}(x|y) = \\dfrac{1}{k} \\sum_{i=1}^{k}(x_i, x_{s_{i,j}})^2
```

- Define the following measure of independence, where ``0 \\leq S \\leq 1``, and 
    low values indicate independence and values close to one occur for 
    synchronized signals.

```math
S^{(k)}(x|y) = \\dfrac{1}{N} \\sum_{i=1}^{N} \\dfrac{R_i^{(k)}(x)}{R_i^{(k)}(x|y)}
```

## Examples

```julia
# Initialise test, specifying embedding dimension, embedding lag 
# and number of nearest neighbors
test = SMeasureTest(m = 4, τ = 3, K = 2:10)
```

## References

1.  Quian Quiroga, R., Arnhold, J. & Grassberger, P. [2000] “Learning 
    driver-response relationships from synchronization patterns,” 
    Phys. Rev. E61(5), 5142–5148.
"""
Base.@kwdef struct SMeasureTest <: CausalityTest
    """ The embedding dimension. Defaults to 2. """
    m::Int = 2
    
    """ The embedding lag. Defaults to 1. """
    τ::Int = 1
    
    """ Neighborhood size(s). Defaults to 2:10. """
    K = 2:10
    
    """ The distance metric used when computing distances. Defaults to `SqEuclidean()`"""
    metric = SqEuclidean()
    
    """ 
    The distance metric used when computing `KDTree` for nearest neighbor searches. 
    Defaults to `Euclidean()`.
    """
    tree_metric = Euclidean()
end

function SMeasureTest(m, τ, K; metric::Metric = SqEuclidean(), tree_metric::Metric = Euclidean())
    SMeasureTest(m = m, τ = τ, K = K; metric = metric, tree_metric = tree_metric)
end

"""

    causality(x, y, test::SMeasureTest)


## Examples 

```julia 
using CausalityTools, DynamicalSystemsBase, Plots

# Create an orbit of the built-in `henon2` map
npts, Ttr = 5000, 500
x, y = columns(trajectory(henon2(c_xy = 1.0), npts, Ttr = Ttr));

# Some random time series
xr, yr = rand(npts), rand(npts);

# Initialise test, specifying embedding dimension, emebdding lag 
# and number of nearest neighbors
Ks = 2:10
test = SMeasureTest(m = 4, τ = 1, K = Ks)

# Compute causality statistic on random time series, and on Henon map
# time series.
Ss_r_xy = causality(xr, yr, test)
Ss_r_yx = causality(yr, xr, test)
Ss_henon_xy = causality(x, y, test)
Ss_henon_yx = causality(y, x, test);

plot(xlabel = "# nearest neighbors (k)", ylabel = "S", ylims = (-0.05, 1.05))
plot!(Ks, Ss_r_xy,  label = "random uncoupled system (x -> y)", marker = stroke(2), c = :black)
plot!(Ks, Ss_r_yx,  label = "random uncoupled system (y -> x)", marker = stroke(2), c = :red)
plot!(Ks, Ss_henon_xy, marker = stroke(2), label = "henon unidir (x -> y)")
plot!(Ks, Ss_henon_yx, marker = stroke(2), label = "henon unidir (y -> x)")
```
"""
function causality(x::AbstractVector, y::AbstractVector, test::SMeasureTest)
    Ss = zeros(Float64, length(test.K))
    
    for (i, k) in enumerate(test.K)
        Ss[i] = s_measure(x, y,test.m, test.τ, k, test.metric, test.tree_metric)
    end
    
    return Ss
end

export SMeasureTest, causality