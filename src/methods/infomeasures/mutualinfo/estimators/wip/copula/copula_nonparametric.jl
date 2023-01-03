import ComplexityMeasures: ProbabilitiesEstimator, DifferentialEntropyEstimator
import StateSpaceSets: AbstractDataset

export Copula
export empirical_copula_transformation

"""
    Copula <: MutualInformationEstimator
    Copula(; est = Kraskov(k = 5), exact = false)

A non-parametric copula-based mutual information estimator.

It is typically many times faster to compute mutual information using `Copula` than
with other [`MutualInformationEstimator`](@ref)s, [`DifferentialEntropyEstimator`](@ref)s,
or [`ProbabilitiesEstimator`](@ref)s, because `Copula` only needs to compute the
entropy of a single (multivariate) variable, whereas the other methods explicitly
computes the entropy of several variables.

If `exact == true`, then the exact empirical cumulative distribution function (ECDF) is
used to compute the empirical copula. If `exact == false`, then a fast sorting-based
approximation to the exact ECDF is computed instead (this breaks ties arbitrarily,
so be careful when applying it to categorical/integer-valued data).

## Description

Assume we have two `Dy`-dimensional and `Dy`-dimensional input [`Dataset`](@ref)s `x` and
`y`, both containing `N` observations. We can define the `Dx + Dy`-dimensional joint
dataset `D = [Dx Dy]`. `Copula` returns the negative *copula entropy* of `D`,
which is equal to the mutual information between `Dx` and `Dy` (Ma & Sun, 2011)[^Ma2011].

[^Ma2011]:
    Ma, J., & Sun, Z. (2011). Mutual information is copula entropy. Tsinghua Science &
    Technology, 16(1), 51-54.

[^Pal2010]:
    Pál, D., Póczos, B., & Szepesvári, C. (2010). Estimation of Rényi entropy and mutual
    information based on generalized nearest-neighbor graphs. Advances in Neural
    Information Processing Systems, 23.
"""
Base.@kwdef struct Copula <: MutualInformationEstimator
    est::Union{ProbabilitiesEstimator, DifferentialEntropyEstimator} = Kraskov(k = 5)
    exact::Bool = false
end

function estimate(measure::MIShannon, est::Copula, x, y)
    X = Dataset(x)
    Y = Dataset(y)
    D = Dataset(X,  Y)
   -entropy(measure.e, est.est, empirical_copula_transformation(D))
end

"""
    rank_transformation(x::AbstractVector)
    rank_transformation(x::AbstractDataset) → ranks::NTuple{D, Vector}

Rank-transform each variable/column of the length-`n` `D`-dimensional dataset `x` and return the
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
function rank_transformation(x::AbstractDataset)
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

"""
    empirical_cdf(x::AbstractVector{<:Real}) → x̄::Vector
    empirical_cdf(x::AbstractDataset) → x̄::Dataset

Rank each sample `xᵢ ∈ x`, rescale it the rank to the interval `[0, 1]` and return
the rescaled ranks `x̄`.
"""
function empcdf(x::AbstractVector)
    F̂ = [count(xᵢ .<= x)  for xᵢ in x] / length(x)
end
empcdf(x::AbstractDataset{D, T}) where {D, T} =
    NTuple{D, Vector{eltype(1.0)}}(empcdf(xⱼ) for xⱼ in columns(x))

# # An example worked out by hand.
# X = Dataset([
#     1 8;
#     2 2;
#     3 6;
#     1 5;
#     2 2;
#     3 1;
#     1 8;
#     2 9;
# ])

# using BenchmarkTools
# using Test
# rank_transform(X)
# @btime rank_transform($X);
# # 2.574 ms (23 allocations: 938.06 KiB)

"""

    empirical_copula_transformation(x::AbstractVector) → empirical_copula::Vector{<:Real}
    empirical_copula_transformation(x::AbstractDataset{D, T}) → empirical_copula::Dataset{D, T}

Apply the empirical copula transformation (as described in Pál et al. (2010)[^Pal2010];
see a summary below) to the each point `xᵢ ∈ x`, where
the `xᵢ` can be either univariate (`x::AbstractVector`) or multivariate
(`x::AbstractDataset`) to compute the empirical copula (here called `empirical_copula)`.

## Description

## Empirical copula transformation

Assume we have a length-`n` sample of data points ``\\bf{X}_{1:n} = \\{\\bf{X}_1, \\bf{X}_2, \\ldots, \\bf{X}_n \\}`` where ``\\bf{X}_i \\in \\mathbb{R}^d``, which is assumed sampled from some distribution ``\\mu`` with density function ``f``. Let ``X_i^j \\in \\mathbb{R}`` denote the j-th coordinate of ``\\bf{X}_i``. Assume these points are represented as the `d`-dimensional [`Dataset`](@ref) which we call `S` (indexed like a matrix where rows are samples).

The *empirical* cumulative distribution function (CDF) for the j-th column of `S`, based on the sample ``\\bf{X}_{1:n}``, is defined as
```math
\\hat{F}_j(y) = \\dfrac{\\left| \\{ 1 \\leq i \\leq n, y \\leq X^j_i \\} \\right|}{n},
```

for any input value ``y \\in \\mathbb{R}`` (which is in general completely unrelated to the j-th column of our sample points). Given the samples ``\\bf{X}_{1:n}``, we can also define a "multivariate empirical CDF" for `S`, ``\\bf{\\hat{F}} : \\mathbb{R}^d \\to [0, 1]^d``, as

```math
\\hat{\\bf{F}}(\\bf{y}) = (\\hat{F}_j(x^1), \\hat{F}_j(x^2), \\ldots, \\hat{F}_j(x^d)),
```

for any point ``\\bf{y} \\in \\mathbb{R}^d`` (which is in general completely unrelated to our sample points, except sharing the property of being `d`-dimensional). Think of this as checking, for each coordinate ``y^j \\in \\bf{y}``, how this coordinate ranks among values in `S[:, j]`.
The map ``\\hat{\\bf{F}}`` is called the *empirical copula transformation*.

Sidenote: Given only our sample, don't actually *know* what the underlying distribution ``\\mu`` is, nor what its cumulative distribution function ``F`` is. But if we did, the analogous map (the *copula transformation*) ``\\bf{F} : \\mathbb{R}^d \\to [0, 1]^d`` would be

```math
\\bf{F}(\\bf{y}) = (F_j(x^1), F_j(x^2), \\ldots, F_j(x^d)).
```

In summary, we've defined the empirical copula *transformation* ``\\hat{\\bf{F}}`` as a map from some `d`-dimensional space to the `d`-dimensional unit square. The j-th axis of ``\\hat{\\bf{F}}``'s domain and the j-th axis of ``\\hat{\\bf{F}}``'s codomain (i.e. the hypersquare) are linked through the *ranks* of `S[:, j]`.

## Empirical copula

The *copula* of ``\\mu`` is the joint distribution ``\\bf{F}(\\bf{X}) = (F_1(X^1), F_2(X^2), \\ldots, F_d(X^d))``. The *empirical copula* (note the lack of "transformation" here) is the set of `d`-dimensional empirical-copula-transformed points ``\\hat{\\bf{Z}} = \\{\\bf{Z}_1, \\bf{Z}_2, \\ldots, \\bf{Z}_n \\} = \\{ \\hat{\\bf{F}}(\\bf{X_1}), \\hat{\\bf{F}}(\\bf{X_2}), \\ldots, \\hat{\\bf{F}}(\\bf{X_n}) \\}``. Note that ``\\hat{\\bf{Z}}`` is an *approximation* of a sample ``\\{\\bf{Z}_1,\\bf{Z}_2, \\ldots, \\bf{Z}_n\\} = \\{\\bf{F}(\\bf{X}_1), \\bf{F}(\\bf{X}_2), \\ldots, \\bf{F}(\\bf{X}_n)\\}`` from the true copula of ``\\mu`` (which we in general don't know, given only some sample points).

[^Pal2010]:
    Pál, D., Póczos, B., & Szepesvári, C. (2010). Estimation of Rényi entropy and mutual
    information based on generalized nearest-neighbor graphs. Advances in Neural
    Information Processing Systems, 23.
"""
function empirical_copula_transformation(x::AbstractDataset{D, T}) where {D, T}
    c = rank_transformation(x) ./ length(x)
    C = Dataset(c...)
end

function empirical_copula_transformation(x::AbstractVector{<:Real})
    rank_transformation(x) ./ length(x)
end

# # An example worked out by hand.
# X = [
#     1 8;
#     2 2;
#     3 6;
#     1 5;
#     2 2;
#     3 1;
#     1 8;
#     2 9;
# ]
# analytical_copula = Dataset([
#     0.125  0.75;
#     0.5    0.25;
#     0.875  0.625;
#     0.25   0.5;
#     0.625  0.375;
#     1.0    0.125;
#     0.375  0.875;
#     0.75   1.0])

# @test copula_transform(Dataset(X)) == analytical_copula
