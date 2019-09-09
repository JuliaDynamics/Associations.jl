
import StatsBase

"""
    VisitationFrequencyTest(k::Int = 1, l::Int = 1, m::Int = 1, n::Int = 1, τ::Int = 1,
        b = 2, 
        estimator::VisitationFrequency = VisitationFrequency(), 
        binning_summary_statistic::Function = StatsBase.mean,
        binning::RectangularBinning, 
        ηs)

The parameters for a transfer entropy test using the `VisitationFrequency` estimator.

## Notes:

- Use `causality(source, target, params::VisitationFrequencyTest)` for regular 
    transfer entropy analysis. This method uses only the `k`, `l`, `m` and ignores `n` 
    when constructing the delay reconstruction. 

- Use `causality(source, target, cond, params::VisitationFrequencyTest)` for conditional 
    transfer entropy analysis. This method uses the `k`, `l`, `m` *and* `n` when constructing 
    the delay reconstruction.

## Example

```
# Prediction lags
ηs = 1:10
binning = RectangularBinning(10)

# Use defaults, binning and prediction lags are required. 
Note that `binning` and `ηs` are *mandatory* keyword arguments.
VisitationFrequencyTest(binning = binning, ηs = ηs)

# The other keywords can also be adjusted
VisitationFrequencyTest(k = 1, l = 2, binning = binning, ηs = ηs)
```
"""
Base.@kwdef struct VisitationFrequencyTest <: TransferEntropyTest
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
    b = 2

    """ 
    If there are several binnings provided, what is the statistic used to summarise the 
    transfer entropy values to a single value?
    """
    binning_summary_statistic::Function = StatsBase.mean

    """ The transfer entropy estimator. """
    estimator::VisitationFrequency = VisitationFrequency()

    """ The binning scheme. """
    binning::Union{RectangularBinning, Vector{RectangularBinning}}

    """ The prediction lags"""
    ηs
end


function causality(source, target, p::VisitationFrequencyTest)
    [p.binning_summary_statistic(
        transferentropy(source, target, p.binning, 
            p.k, p.l, p.m, η = η, τ = p.τ, 
            estimator = p.estimator)) for η in p.ηs]
end

export VisitationFrequencyTest