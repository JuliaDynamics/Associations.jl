using Neighborhood: searchstructure, bulkisearch
using Neighborhood: NeighborNumber, Theiler
using StateSpaceSets: AbstractStateSpaceSet
using Distances: SqEuclidean, Euclidean
using Distances: pairwise, evaluate

export SMeasure
export s_measure

"""
    SMeasure < AssociationMeasure
    SMeasure(; K::Int = 2, dx = 2, dy = 2, τx = - 1, τy = -1, w = 0)

`SMeasure` is a bivariate association measure from [Arnhold1999](@citet)
and [Quiroga2000](@citet) that measure directional dependence
between two input (potentially multivariate) time series.

Note that `τx` and `τy` are negative; see explanation below.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for directional dependence.
- Use with [`s_measure`](@ref) to compute the raw s-measure statistic.

## Description

The steps of the algorithm are:

0. From input time series ``x(t)`` and ``y(t)``, construct the delay embeddings (note
    the positive sign in the embedding lags; therefore inputs parameters
    `τx` and `τy` are by convention negative).

```math
\\begin{align*}
\\{\\bf{x}_i \\} &= \\{(x_i, x_{i+\\tau_x}, \\ldots, x_{i+(d_x - 1)\\tau_x}) \\} \\\\
\\{\\bf{y}_i \\} &= \\{(y_i, y_{i+\\tau_y}, \\ldots, y_{i+(d_y - 1)\\tau_y}) \\} \\\\
\\end{align*}
```

1. Let ``r_{i,j}`` and ``s_{i,j}`` be the indices of the `K`-th nearest neighbors
    of ``\\bf{x}_i `` and ``\\bf{y}_i``, respectively. Neighbors closed than `w` time indices
    are excluded during searches (i.e. `w` is the Theiler window).

2. Compute the the mean squared Euclidean distance to the ``K`` nearest neighbors
    for each ``x_i``, using the indices ``r_{i, j}``.

```math
R_i^{(k)}(x) = \\dfrac{1}{k} \\sum_{i=1}^{k}(\\bf{x}_i, \\bf{x}_{r_{i,j}})^2
```

- Compute the y-conditioned mean squared Euclidean distance to the ``K`` nearest
    neighbors for each ``x_i``, now using the indices ``s_{i,j}``.

```math
R_i^{(k)}(x|y) = \\dfrac{1}{k} \\sum_{i=1}^{k}(\\bf{x}_i, \\bf{x}_{s_{i,j}})^2
```

- Define the following measure of independence, where ``0 \\leq S \\leq 1``, and
    low values indicate independence and values close to one occur for
    synchronized signals.

```math
S^{(k)}(x|y) = \\dfrac{1}{N} \\sum_{i=1}^{N} \\dfrac{R_i^{(k)}(x)}{R_i^{(k)}(x|y)}
```

## Input data

The algorithm is slightly modified from [Grassberger1999](@cite) to allow univariate timeseries as input.

- If `x` and `y` are [`StateSpaceSet`](@ref)s then use `x` and `y` as is and ignore the parameters
    `dx`/`τx` and `dy`/`τy`.
- If `x` and `y` are scalar time series, then create `dx` and `dy` dimensional embeddings,
    respectively, of both `x` and `y`, resulting in `N` different `m`-dimensional embedding points
    ``X = \\{x_1, x_2, \\ldots, x_N \\}`` and ``Y = \\{y_1, y_2, \\ldots, y_N \\}``.
    `τx` and `τy` control the embedding lags for `x` and `y`.
- If `x` is a scalar-valued vector and `y` is a [`StateSpaceSet`](@ref), or vice versa,
    then create an embedding of the scalar timeseries using parameters `dx`/`τx` or `dy`/`τy`.

In all three cases, input StateSpaceSets are length-matched by eliminating points at the end of
the longest StateSpaceSet (after the embedding step, if relevant) before analysis.
"""
Base.@kwdef struct SMeasure{M, TM} <: AssociationMeasure
    K::Int = 2
    metric::M = SqEuclidean()
    tree_metric::TM = Euclidean()
    τx::Int = -1
    τy::Int = -1
    dx::Int = 2
    dy::Int = 2
    w::Int = 0
end

"""
    s_measure(measure::SMeasure, x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet)

Compute the [`SMeasure`](@ref) from source `x` to target `y`.
"""
function s_measure(measure::SMeasure, x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet)
    return estimate(measure, x, y)
end

# Internal method for use with `independence`
function estimate(measure::SMeasure, x::AbstractStateSpaceSet, y::AbstractStateSpaceSet)
    (; K, metric, tree_metric, τx, τy, dx, dy, w) = measure

    # Match length of StateSpaceSets by excluding end points.
    lx = length(x); ly = length(y)
    lx > ly ? X = x[1:ly, :] : X = x
    ly > lx ? Y = y[1:lx, :] : Y = y
    N = length(X)

    T = eltype(1.0)
     # Pre-allocate vectors to hold indices and distances during loops
    dists_x = zeros(T, K)
    dists_x_cond_y = zeros(T, K)

    # Mean squared distances in X, and
    # mean squared distances in X conditioned on Y
    Rᵢᵏx = zeros(T, N)
    Rᵢᵏxy = zeros(T, N)

    # Search for the K nearest neighbors of each points in both X and Y
    treeX = searchstructure(KDTree, X, tree_metric)
    treeY = searchstructure(KDTree, Y, tree_metric)
    neighborhoodtype, theiler = NeighborNumber(K), Theiler(w)
    idxs_X = bulkisearch(treeX, X, neighborhoodtype, theiler)
    idxs_Y = bulkisearch(treeY, Y, neighborhoodtype, theiler)

    for n in 1:N
        pxₙ = X[n]

        for j = 1:K
            rₙⱼ = idxs_X[n][j] # nearest neighbor indices in X
            sₙⱼ = idxs_Y[n][j] # nearest neighbor indices in Y
            dists_x[j] = evaluate(metric, pxₙ, X[rₙⱼ])
            dists_x_cond_y[j] = evaluate(metric, pxₙ, X[sₙⱼ])
        end

        Rᵢᵏx[n] = sum(dists_x) / K
        Rᵢᵏxy[n] = sum(dists_x_cond_y) / K
    end

    return sum(Rᵢᵏx ./ Rᵢᵏxy) / N
end
