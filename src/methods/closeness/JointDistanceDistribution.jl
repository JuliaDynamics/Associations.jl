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
    JointDistanceDistribution(; metric = Euclidean(), B = 10, D = 2, τ = 1, μ = 0.0)

The joint distance distribution (JDD) measure (Amigó & Hirata, 2018)[^Amigo2018].

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for directional dependence.
- Use with [`jdd`](@ref) to compute the joint distance distribution `Δ` from Amigó & Hirata (2018).

## Keyword arguments

- **`distance_metric::Metric`**: An instance of a valid distance metric from `Distances.jl`.
    Defaults to `Euclidean()`.
- **`B::Int`**: The number of equidistant subintervals to divide the interval `[0, 1]` into
    when comparing the normalised distances.
- **`D::Int`**: Embedding dimension.
- **`τ::Int`**: Embedding delay.
- **`μ`**: The hypothetical mean value of the joint distance distribution if there
    is no coupling between `x` and `y` (default is `μ = 0.0`).


## Examples

* [Computing the JJD](@ref quickstart_jdd)
* [Independence testing using JJD](@ref quickstart_jddtest)

[^Amigo2018]:
    Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate
    flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of
    Nonlinear Science 28.7 (2018): 075302.
"""
Base.@kwdef struct JointDistanceDistribution{M, T} <: AssociationMeasure
    metric::M = Euclidean()
    B::Int = 5
    D::Int = 3
    τ::Int = 1
    μ::T = 0.0
end

# The convenience wrapper `jdd`` is in deprecations folder for now.

function estimate(measure::JointDistanceDistribution, est::Nothing, source, target)
    return estimate(measure, source, target)
end

# Internal method for compatibility with independence tests.
function estimate(measure::JointDistanceDistribution, source, target)
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
