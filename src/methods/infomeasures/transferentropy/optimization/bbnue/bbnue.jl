using DelayEmbeddings, Statistics, TimeseriesSurrogates
include("candidate_variable_sets.jl")

export bbnue

# Keep this for when a common interface for optimized variable selection methods has been established
# """
#export BBNUE
#     BBNUE(est) <: TransferDifferentialInfoEstimator

# The bootstrap-based non-uniform embedding estimator (BB-NUE) for conditional transfer entropy (Montalto et al., 2014).
# Uses the estimator `est` to compute relevant marginal entropies. (e.g. `VisitationFrequency(RectangularBinning(3))`)

# [^Montalto2014]: Montalto, A.; Faes, L.; Marinazzo, D. MuTE: A MATLAB toolbox to compare established and novel estimators of the multivariate transfer entropy. PLoS ONE 2014, 9, e109462.
# """
# struct BBNUE{E} <: TransferDifferentialInfoEstimator
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
which is what Montalto et al. (2014)[^Montalto2014] uses for their bin-estimator.
This implementation accepts *any* [`DifferentialInfoEstimator`](@ref) or
[`ProbabilitiesEstimator`](@ref) that accepts multivariate data or that implements
[`marginal_encodings`](@ref).

Montalto et al.'s bin-estimator corresponds to using the `VisitationFrequency` estimator
with bins whose sides are equal-length, e.g. `VisitationFrequency(RectangularBinning(0.5))`.
In this implementation, any rectangular binning can be used.

## Input data

Multivariate `source`, `target` and `cond` (if given) can be given as univariate
`AbstractVector`s or as multivariate `StateSpaceSet`s or `Vector{AbstractVector}`.
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
bbnue(est, s, t, StateSpaceSet(c1, c2, c3))
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
using CausalityTools
using DynamicalSystemsBase: trajectory
sys = logistic2_unidir(c_xy = 0.8, r₁ = 3.78, r₂ = 3.92) # built-in system
n = 5000
orbit = trajectory(sys, n, Ttr = 10000)
x, y = columns(orbit)
z, w = randn(n), rand(n) # some completely unrelated timeseries to condition on.

# Use a coarse-grained rectangular binning with subdivisions in each dimension,
# to keep computation costs low and to ensure the probability distributions
# over the bins don't approach the uniform distribution (need enough points
# to fill bins).
est = ValueHistogram(4)

te_xy = bbnue(est, x, y, surr = RandomShuffle(), nsurr = 50)
te_yx = bbnue(est, y, x, surr = RandomShuffle(), nsurr = 50)
te_xy, te_yx
```

[^Montalto2014]:
    Montalto, A.; Faes, L.; Marinazzo, D. MuTE: A MATLAB toolbox to compare established
    and novel estimators of the multivariate transfer entropy. PLoS ONE 2014, 9, e109462.
"""
function bbnue(measure::TransferEntropy, est::DiscreteOrDifferentialInfoEstimator, x...;
        η = 1, include_instantaneous = true,
        method_delay = "ac_min", maxlag::Union{Int, Float64} = 0.05,
        α = 0.05, nsurr = 19, surr::TimeseriesSurrogates.Surrogate = RandomShuffle())

    processed_inputs = process_input.(x)
    Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond =
        candidate_embedding(processed_inputs...;
            η,
            include_instantaneous,
            method_delay,
            maxlag)
    return optim_te(measure, est, Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond;
        α, nsurr, surr)
end

bbnue(est::DiscreteOrDifferentialInfoEstimator, x...; kwargs...) =
    bbnue(TEShannon(), est, x...; kwargs...)


function optim_te(measure::TransferEntropy, est::DiscreteOrDifferentialInfoEstimator,
        Ω, Y⁺, τs, js, idxs_source, idxs_target, idxs_cond;
        α = 0.05, nsurr = 100, surr::Surrogate = RandomShuffle())
    e = measure.e
    τs_comb = [(τs...)...,]
    js_comb = [(js...)...,]

    n_candidate_variables = length(Ω)

    𝒮 = Vector{Vector{Float64}}(undef, 0)
    𝒮_τs = Vector{Int}(undef, 0)
    𝒮_js = Vector{Int}(undef, 0)

    k = 1
    while k <= n_candidate_variables
        n_remaining_candidates = length(Ω)
        CMIs_between_Y⁺_and_candidates = zeros(n_remaining_candidates)

        # At first iteration, only loop through source variable. If no source variable is found that
        # yields significant TE, terminate.
        for i = 1:n_remaining_candidates
            if k == 1 || length(𝒮) == 0
                Cᵢ = Ω[i]
                CMI_Y⁺_Cᵢ =
                    entropy(e, est, StateSpaceSet(Y⁺, StateSpaceSet(Cᵢ))) -
                    entropy(e, est, StateSpaceSet(Cᵢ))
            else
                Cᵢ = [Ω[i], 𝒮...]
                CMI_Y⁺_Cᵢ =
                    entropy(e, est, StateSpaceSet(Y⁺, StateSpaceSet(Cᵢ...,))) -
                    entropy(e, est, StateSpaceSet(Cᵢ...,))
            end
            CMIs_between_Y⁺_and_candidates[i] = CMI_Y⁺_Cᵢ
        end

        idx = findfirst(x -> x == minimum(CMIs_between_Y⁺_and_candidates), CMIs_between_Y⁺_and_candidates)
        cₖ = Ω[idx]

        # Test the significance of this candidate by using a permutation test. The type of surrogate
        # is given by `surr`, and we will use `nsurr` surrogate realizations.
        CMI_permutations = zeros(nsurr)
        s = surrogenerator(cₖ, surr)

        # If k == 1, no candidates have been selected, so CMI reduces to MI
        if k == 1
            condmutualinfoₖ = CMIs_between_Y⁺_and_candidates[idx]

            for i = 1:nsurr
                surr_cₖ = s() # Surrogate version of cₖ
                CMI_permutations[i] = mutualinfo(est, Y⁺, surr_cₖ)
            end
        # If k > 1, at least one candidate has been selected, so we compute CMI
        else
            # Precompute terms that do not change during surrogate loop
            H_Y⁺_𝒮 = entropy(e, est, StateSpaceSet(Y⁺, StateSpaceSet(𝒮...,)))
            H_𝒮 = entropy(e, est, StateSpaceSet(𝒮...))

            # Original TE
            condmutualinfoₖ = H_Y⁺_𝒮 +
                    entropy(e, est, StateSpaceSet([cₖ, 𝒮...,]...,)) -
                    entropy(e, est, StateSpaceSet(Y⁺, StateSpaceSet([cₖ, 𝒮...,]...,))) -
                    H_𝒮

            for i = 1:nsurr
                surr_cₖ = s() # Surrogate version of cₖ
                CMI_permutations[i] = H_Y⁺_𝒮 +
                    entropy(e, est, StateSpaceSet([surr_cₖ, 𝒮...]...,)) -
                    entropy(e, est, StateSpaceSet(Y⁺, StateSpaceSet([surr_cₖ, 𝒮...]...,))) -
                    H_𝒮
            end
        end

        # If the candidate passes the significance test, add it to list of selected candidates
        # and remove it from list of remaining candidates.
        if condmutualinfoₖ > quantile(CMI_permutations, 1 - α)
            push!(𝒮, cₖ)
            push!(𝒮_τs, τs_comb[idx])
            push!(𝒮_js, js_comb[idx])
            deleteat!(Ω, idx)
            deleteat!(τs_comb, idx)
            deleteat!(js_comb, idx)
            k = k + 1
        # If the candidate does not pass significance test, terminate.
        else
            k = n_candidate_variables + 1
        end
    end

    # No variables were selected at all.
    if length(𝒮) == 0
        return 0.0, Int[], Int[], idxs_source, idxs_target, idxs_cond
    end

    # If no variables were selected from the source process, then TE is not well-defined.
    n_source_vars_picked = count(x -> x ∈ idxs_source, 𝒮_js)
    if n_source_vars_picked == 0
        return 0.0, Int[], Int[], idxs_source, idxs_target, idxs_cond
    end

    # No variables were selected from the target or conditional processes.
    𝒮_nonX = [ts for (ts, j) in zip(𝒮, 𝒮_js) if j ∉ idxs_source]
    if length(𝒮_nonX) == 0
        return 0.0, Int[], Int[], idxs_source, idxs_target, idxs_cond
    end

    CE2 = entropy(e, est, StateSpaceSet(Y⁺, StateSpaceSet(𝒮...,))) -
        entropy(e, est, StateSpaceSet(𝒮...,))

    CE1 = entropy(e, est, StateSpaceSet(Y⁺, StateSpaceSet(𝒮_nonX...,))) -
        entropy(e, est, StateSpaceSet(𝒮_nonX...,))

    CMI = CE1 - CE2
    return CMI, 𝒮_js, 𝒮_τs, idxs_source, idxs_target, idxs_cond
end
