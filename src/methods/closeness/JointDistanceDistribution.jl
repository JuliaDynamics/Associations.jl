using HypothesisTests
using Distances: Euclidean
using Distances: pairwise
using DelayEmbeddings: genembed

import HypothesisTests: OneSampleTTest, pvalue
export OneSampleTTest, pvalue
export JointDistanceDistribution
export jdd

function normalise_minmax(x::T, vmin, vmax) where T
    if x == zero(T)
        return zero(T)
    else
        return (x - vmin)/(vmax - vmin)
    end
end

"""
    JointDistanceDistribution <: AssociationMeasure end
    JointDistanceDistribution(; metric = Euclidean(), B = 10, D = 2, τ = -1, μ = 0.0)

The joint distance distribution (JDD) measure [Amigo2018](@cite).

## Usage

- Use with [`association`](@ref)/[`jdd`](@ref) to compute the joint distance distribution measure `Δ` from
    [Amigo2018](@citet).
- Use with [`independence`](@ref) to perform a formal hypothesis test for directional
    dependence.

## Keyword arguments

- **`distance_metric::Metric`**: An instance of a valid distance metric from `Distances.jl`.
    Defaults to `Euclidean()`.
- **`B::Int`**: The number of equidistant subintervals to divide the interval `[0, 1]` into
    when comparing the normalised distances.
- **`D::Int`**: Embedding dimension.
- **`τ::Int`**: Embedding delay. By convention, `τ` is negative.
- **`μ`**: The hypothetical mean value of the joint distance distribution if there
    is no coupling between `x` and `y` (default is `μ = 0.0`).

## Description

From input time series ``x(t)`` and ``y(t)``, we first construct the delay embeddings (note
the positive sign in the embedding lags; therefore the input parameter
`τ` is by convention negative).

```math
\\begin{align*}
\\{\\bf{x}_i \\} &= \\{(x_i, x_{i+\\tau}, \\ldots, x_{i+(d_x - 1)\\tau}) \\} \\\\
\\{\\bf{y}_i \\} &= \\{(y_i, y_{i+\\tau}, \\ldots, y_{i+(d_y - 1)\\tau}) \\} \\\\
\\end{align*}
```

The algorithm then proceeds to analyze the distribution of distances between points
of these embeddings, as described in [Amigo2018](@citet).

## Examples

* [Computing the JDD](@ref quickstart_jdd)
* [Independence testing using JDD](@ref quickstart_jddtest)
"""
Base.@kwdef struct JointDistanceDistribution{M, T} <: AssociationMeasure
    metric::M = Euclidean()
    B::Int = 5
    D::Int = 3
    τ::Int = -1
    μ::T = 0.0
end

# The convenience wrapper `jdd`` is in deprecations folder for now.

function association(measure::JointDistanceDistribution, est::Nothing, source, target)
    return association(measure, source, target)
end

# Internal method for compatibility with independence tests.
function association(measure::JointDistanceDistribution, source, target)
    (; metric, B, D, τ) = measure
    length(source) == length(target) || error("lengths of inputs must match")
    js = ([1 for i = 1:D]...,)
    τs = (collect(0:-τ:-(D-1)*τ)...,)
    Ex = genembed(source, τs, js)
    Ey = genembed(target, τs, js)
    Dx = pairwise(metric, Ex.data)
    Dy = pairwise(metric, Ey.data)

    # Normalise the distances to the interval [0, 1]
    fDx = filter(dx -> dx .> 0, Dx)
    fDy = filter(dy -> dy .> 0, Dy)
    Dx_min, Dx_max = minimum(fDx), maximum(fDx)
    Dy_min, Dy_max = minimum(fDy), maximum(fDy)

    # We don't simply normalize to [0, 1]. We transform the distances to a uniform distribution
    # over [0, 1] using the normalized rank transformation.
    normDx = rank_transformation(vec(Dx)) ./ length(Dx)
    normDy = rank_transformation(vec(Dy)) ./ length(Dy)

    δs = fill(1.0, 2 * B)
    for (k, b) in enumerate(1:(2 * B))
        bmin = (b - 1) / 2B
        bmax = b / 2B

        δ̃min = jdd_step3(normDx, normDy, bmin, bmax)
        δs[k] = δ̃min
    end

    Δ = [δs[B + i] - δs[i] for i in 1:B]
    return Δ
end

function jdd_step3(Dx, Dy, bmin, bmax)
    δy_min = 1.0 # by definition, 1.0 is the largest normalized distance.
    found = false
    for (i, δx) in enumerate(Dx)
        if bmin <= δx < bmax
            δy = Dy[i]
            if δy < δy_min
                δy_min = δy
                found = true
            end
        end
    end
    if found
        return δy_min
    else
        return bmin
    end
end

"""
    rank_transformation(x::AbstractVector)
    rank_transformation(x::AbstractStateSpaceSet) → ranks::NTuple{D, Vector}

Rank-transform each variable/column of the length-`n` `D`-dimensional StateSpaceSet `x` and return the
rank-transformed variables as a `D`-tuple of length-`n` vectors.

Returns the unscaled `ranks`. Divide by `n` to get an *approximation* to the
empirical cumulative distribution function (ECDF)  `x`.

## Description

Modulo division by `n`, `rank_transformation` does *roughly* the same as naively computing the ECDF as
```julia
[count(xᵢ .<= x)  for xᵢ in x] / length(x)
```

but an order of magnitude faster and with roughly three orders of magnitude less
allocations. The increased efficiency of this function relative to naively computing the
ECDF is
because it uses sorting of the input data to determine ranks,
arbitrarily breaking ties according to the sorting algorithm. Rank ties can therefore
never occur, and equal values are assigned different but close ranks. To preserve
ties, which you might want to do for example when dealing with
categorical or integer-valued data, use (the much slower) [`empcdf`](@ref).
"""
function rank_transformation(x::AbstractStateSpaceSet)
    s = zeros(Int, length(x)) # re-use for each marginal
    [rank_transformation!(s, xⱼ) for xⱼ in columns(x)]
end

function rank_transformation(x::AbstractVector{T}) where T
    N = length(x)
    s = zeros(Int, N)
    return rank_transformation!(s, x)
end

function rank_transformation!(
        s::AbstractVector{Int},
        x::AbstractVector{T}) where T <: Real
    N = length(x)
    r = zeros(N)
    # Break ties arbitrarily by sorting. This means that ties are broken according to the
    # sorting algorithm used, and equal values are assigned different ranks.
    sortperm!(s, x)
    for j in 1:N
        r[s[j]] = j
    end
    return r
end