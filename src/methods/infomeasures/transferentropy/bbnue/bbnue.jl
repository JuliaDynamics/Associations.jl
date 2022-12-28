using TimeseriesSurrogates
include("autoutils.jl")

export bbnue

# Keep this for when a common interface for optimized variable selection methods has been established
# """
#export BBNUE
#     BBNUE(est) <: TransferDifferentialEntropyEstimator

# The bootstrap-based non-uniform embedding estimator (BB-NUE) for conditional transfer entropy (Montalto et al., 2014).
# Uses the estimator `est` to compute relevant marginal entropies. (e.g. `VisitationFrequency(RectangularBinning(3))`)

# [^Montalto2014]: Montalto, A.; Faes, L.; Marinazzo, D. MuTE: A MATLAB toolbox to compare established and novel estimators of the multivariate transfer entropy. PLoS ONE 2014, 9, e109462.
# """
# struct BBNUE{E} <: TransferDifferentialEntropyEstimator
#     est::E
# end

"""
    bbnue(source, target, [cond], est;
        η = 1, include_instantaneous = true,
        method_delay = "ac_min", maxlag::Union{Int, Float64} = 0.05,
        surr::Surrogate = RandomShuffle(), nsurr = 100, α = 0.05)
        ) → te, js, τs, idxs_source, idxs_target, idxs_cond

Estimate transfer entropy (TE) from `source` to `target` (conditioned on `cond` if given)
for prediction lag `η`, using the bootstrap-based non-uniform embedding (BBNUE)
method from Montalta et al. (2004) [^Montalto2014].

## Implementation details

The BBNUE method uses a bootstrap-based criterion to identify the most relevant and
minimally redundant variables from the the past of `target`, present/past of `source`,
and (if given) the present/past of `cond`, that contribute most to `target`'s future.

This implementation uses a conditional entropy minimization criterion for selecting variables,
which is what Montalto et al. (2014)[^Montalto2014] uses for their bin-estimator. Here, any
estimator can be used, but variables will be selected using conditional entropy minimization
regardless of the choice of estimator.

Montalto et al.'s bin-estimator corresponds to using the `VisitationFrequency` estimator with bins
whose sides are equal-length, e.g. `VisitationFrequency(RectangularBinning(0.5))`.
In this implementation, any rectangular binning can be used.

## Input data

Multivariate `source`, `target` and `cond` (if given) can be given as univariate
`AbstractVector`s or as multivariate `Dataset`s or `Vector{AbstractVector}`.

For example, if you want to compute
the BBNUE-transfer entropy from a univariate source to a univariate target,
potentially conditioned on many different variables, you can do the following:

```julia
n = 1000
# Source and target variables
s, t = rand(n), rand(n)

# Variables that might potentially influence `t` along with `s`
c1, c2, c3 = rand(n), rand(n), rand(n)

est = NaiveKernel(0.3)
bbnue(s, t, Dataset([c1, c2, c3]), est)
```

## Variable selection and significance testing

In this implementation, the maximum lag for each embedding variable is determined using `estimate_delay`
from `DelayEmbeddings`. The keywords `method_delay` (the default is `"ac_min"`) controls the method
for estimating the delay, and `maxlag` is the maximum allowed delay (if `maxlag ∈ [0, 1]` is a fraction,
then the maximum lag is that fraction of the input time series length, and if `maxlag` is an integer,
then the maximum lag is `maxlag`).

If `instantaneous` is `true`, then instantaneous interactions are also considered, i.e. effects like
`source(t) → target(t)` are allowed. `η` is the forward prediction lag.

Significance testing is performed using a permutation test. At each iteration of the
variable selection routine, we first compute the transfer entropy using the new candidate
``c_k``. Then, the computation is repeated `nsurr` times, at each iteration replacing ``c_k``
with a surrogate of type `surr`. If transfer entropy using the original ``c_k`` exceeds the
the `1 - α`-quantile of that of the surrogate ensemble, then ``c_k`` is deemed significant
to the future of `target` and is included in the set of selected variables.

If no relevant variables pass the permutation test, then TE is not well-defined, and a value of `0.0`
is returned.

## Returns

A 6-tuple is returned, consisting of:
- `te`: The computed transfer entropy value.
- `js`: The indices of the selected variables. `js[i]` is the `i`-th entry in the array `[idxs_source..., idxs_target..., idxs_cond...,]`.
- `τs`: The embedding lags of the selected variables. `τs[i]` corresponds to `js[i]`.
- `idxs_source`: The indices of the source variables.
- `idxs_target`: The indices of the target variables.
- `idxs_cond`: The indices of the conditional variables (empty if `cond` is not given).

## Example

```julia
using CausalityTools, DynamicalSystems
sys = ExampleSystems.logistic2_unidir(c_xy = 0.8, r₁ = 3.78, r₂ = 3.92)
orbit = trajectory(sys, 10000, Ttr = 10000)
x, y = columns(orbit)

# Use a coarse-grained rectangular binning with subdivisions in each dimension,
# to keep computation costs low and to ensure the probability distributions
# over the bins don't approach the uniform distribution (need enough points
# to fill bins).
est = NaiveKernel(0.3)
te_xy = bbnue(x, y, est, surr = RandomShuffle(), nsurr = 100, include_instantaneous = true)
te_yx = bbnue(y, x, est, surr = RandomShuffle(), nsurr = 100, include_instantaneous = true)

te_xy, te_yx
```

[^Montalto2014]: Montalto, A.; Faes, L.; Marinazzo, D. MuTE: A MATLAB toolbox to compare established and novel estimators of the multivariate transfer entropy. PLoS ONE 2014, 9, e109462.
"""
function bbnue(e::EntropyDefinition, est::ProbabilitiesEstimator, source, target, cond;
        η = 1, include_instantaneous = true,
        method_delay = "ac_min", maxlag::Union{Int, Float64} = 0.05,
        α = 0.05, nsurr = 19, surr::TimeseriesSurrogates.Surrogate = RandomShuffle())

    Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond =
        candidate_embedding(
            process_input(source),
            process_input(target),
            process_input(cond);
            η,
            include_instantaneous,
            method_delay,
            maxlag)

    return optim_te(e, Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond, est;
        α, nsurr, surr)
end

function bbnue(e::EntropyDefinition, est::ProbabilitiesEstimator, source, target;
        η = 1, include_instantaneous = true,
        method_delay = "ac_min", maxlag::Union{Int, Float64} = 0.05,
        α = 0.05, nsurr = 19, surr::TimeseriesSurrogates.Surrogate = RandomShuffle())

    Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond =
        candidate_embedding(
            process_input(source),
            process_input(target);
            η,
            include_instantaneous,
            method_delay,
            maxlag)

    return optim_te(e, est, Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond;
        α, nsurr, surr)
end

bbnue(est::ProbabilitiesEstimator, s, t; base = 2, kwargs...) =
    bbnue(Shannon(; base), est, s, t; kwargs...)

bbnue(est::ProbabilitiesEstimator, s, t, c; base = 2, kwargs...) =
    bbnue(Shannon(; base), est, s, t,c ; kwargs...)
