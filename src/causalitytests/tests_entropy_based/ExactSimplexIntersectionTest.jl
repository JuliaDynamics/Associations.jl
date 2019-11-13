import StatsBase
import TransferEntropy: TEVars, transferentropy, te_embed
import PerronFrobenius: AbstractTriangulationInvariantMeasure, invariantmeasure
import CausalityToolsBase: TriangulationBinning, ExactIntersection

"""
    ExactSimplexIntersectionTest

The parameters for a transfer entropy test where the invariant measure is estimated
using an approximation to the transfer operator over a triangulation of the delay
reconstruction [1]. Essentially, this allows the generation of an unlimited amount
of points from the transfer operator. Transfer entropy then can be estimated by
superimposing a rectangular grid as usual over the points generated from the 
transfer operator.

## Notes

- Compared to the rectangular estimators, the exact simplex intersection approach is much slower.
    However, it might be beneficial for sparse time series[1].
    A faster triangulation estimator is the `ApproximateSimplexIntersectionTest`.
- If you're doing any sort of sensitivity analysis over binning schemes, then use the low-level
    estimators to first construct the invariant measure over the triangulation. You can then
    estimate transfer entropy "for free" from the points you generate from the triangulation.
    Over these points you can superimpose any grid you want and quickly compute transfer entropy
    using one of the rectangular estimators.


## Mandatory keyword arguments

- **`binning::RectangularBinning`**: An instance of a [`RectangularBinning`](@ref) that dictates
    how the delay embedding is discretized.

- **`ηs`**: The prediction lags (that gos into the ``T_{f}`` component of the embedding).

## Optional keyword arguments

- **`k::Int`**: The dimension of the ``T_{f}`` component of the embedding.

- **`l::Int`**: The dimension of the ``T_{pp}`` component of the embedding.

- **`m::Int`**: The dimension of the ``S_{pp}`` component of the embedding.

- **`n::Int`**: The dimension of the ``C_{pp}`` component of the embedding.

- **`τ::Int`**: The embedding lag. Default is `τ = 1`.

- **`b`**: Base of the logarithm. The default (`b = 2`) gives the TE in bits.

- **`estimator::Union{VisitationFrequency, TransferOperatorGrid}`**: A `VisitationFrequency`
    or a `TransferOperatorGrid` estimator instance.

- **`n_pts::Int`**: The number of points to generate from the invariant distribution
    over the triangulation. It is over these points transfer entropy is finally
    generated using the provided rectangular `estimator`. Defaults to `n_pts = 10000`.

- **`binning_summary_statistic::Function`**: A summary statistic to summarise the
    transfer entropy values if multiple binnings are provided.

## Examples

```julia
using DynamicalSystems

# Short time series of two coupled logistic maps (x drives y)
sys = logistic2_unidir(c_xy = 0.5)
npts = 75
x, y = columns(trajectory(sys, npts, Ttr = 1000))

# average over a few different binning schemes, use a few different
# prediction lags and use the transfer operator grid estimator on the point
# cloud generated from the invariant measure over the triangulation.
binnings = [RectangularBinning(i) for i = 3:5]
estimator = TransferOperatorGrid()
ηs = 1:3

# Perform causality test, both from x to y and y to x
test = ExactSimplexIntersectionTest(binning = binnings, estimator = estimator, ηs = ηs)
te_xtoy = causality(x, y, test)
te_ytox = causality(y, x, test)
```

## References

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
    [https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.042212](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.042212)
"""
Base.@kwdef struct ExactSimplexIntersectionTest{N} <: TransferEntropyCausalityTest{N}
    """ The delay reconstruction parameter k (controls dimension of ``T_{f}`` component of embedding). """
    k::Int = 1

    """ The delay reconstruction parameter l (controls dimension of ``T_{pp}`` component of embedding). """
    l::Int = 1

    """ The delay reconstruction parameter m (controls dimension of ``S_{pp}`` component of embedding). """
    m::Int = 1

    """ The delay reconstruction parameter n (controls dimension of ``C_{pp}`` component of embedding). """
    n::Int = 1

    """ The delay reconstruction lag for the ``T_{pp}`` component of the embedding. """
    τ::Int = 1

    """ The base of the logarithm for computing TE. """
    b::Number = 2

    """
    The transfer entropy estimator used to estimate transfer entropy *after* the invariant measure
    over the triangulated delay reconstruction has been estimated.
    """
    estimator::TransferEntropyEstimator = TransferOperatorGrid()

    """ The number of points to generate from the invariant distribution over the triangulation. """
    n_pts::Int = 10000

    """
    If there are several binnings provided, what is the statistic used to summarise the
    transfer entropy values to a single value?
    """
    binning_summary_statistic::Function = StatsBase.mean

    """
    The binning scheme(s). If more than one is provided, the `binning_summary_statistic` is
    applied to the computed transfer entropy values, and a single value is returned.
    """
    binning::Union{RectangularBinning, Vector{RectangularBinning}}


    """ The prediction lags"""
    ηs

    function ExactSimplexIntersectionTest(k::Int, l::Int, m::Int, n::Int, τ::Int, b::Number, 
            estimator::E, n_pts::Int, 
            binning_summary_statistic::Function, 
            binning::Union{RectangularBinning, Vector{RectangularBinning}}, 
            ηs) where {E <: TransferEntropyEstimator}

        N = length(ηs) # length of return vector when used with `causality`
        return new{N}(k, l, m, n, τ, b, estimator, n_pts, 
            binning_summary_statistic, binning, ηs)
    end
end


function causality(source::AbstractVector{T}, target::AbstractVector{T},
        test::ExactSimplexIntersectionTest)  where {T <: Real}

    k, l, m, τ, ηs = test.k, test.l, test.m, test.τ, test.ηs
    k + l + m >= 3 || throw(ArgumentError("`dim = k + l + m` must be 3 or higher for regular TE"))

    # Delay reconstructions and constructing TEVars instances, one for each η,
    # then compute invariant measure over triangulation once for each η
    te_vars = Vector{TEVars}(undef, length(test.ηs))
    invariant_measures_over_triang = Vector{AbstractTriangulationInvariantMeasure}(undef, length(test.ηs))

    for (i, η) in enumerate(ηs)
        pts, vars = te_embed(source, target, k, l, m, η = η, τ = τ)
        te_vars[i] = vars

        # ugly and slow hack, PerronFrobenius needs to accept SVectors and Datasets by default.
        regpoints = [Vector(x) for x in pts.reconstructed_pts]

        # Compute invariant measure over triangulation
        mu = invariantmeasure(regpoints, TriangulationBinning(), ExactIntersection())
        invariant_measures_over_triang[i] = mu
    end


    # Superimpose rectangular grid on the triangulations and compute transfer entropy
    # from that
    tes_over_ηs = zeros(Float64, length(ηs))

    for i = 1:length(ηs)
        mu = invariant_measures_over_triang[i]
        tes_over_binnings = [transferentropy(mu, te_vars[i], binning, estimator = test.estimator, n = test.n_pts, b = test.b)
                for binning in test.binning]
        tes_over_ηs[i] = test.binning_summary_statistic(tes_over_binnings)
    end

    return tes_over_ηs
end

export ExactSimplexIntersectionTest
