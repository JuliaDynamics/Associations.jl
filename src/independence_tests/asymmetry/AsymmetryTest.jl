using ShiftedArrays
export AsymmetryTest

include("dispatch.jl")

"""
    AsymmetryTest <: IndependenceTest
    AsymmetryTest(measure, [est];
        ηT = 1:10
        τS = 0
        τT = 0:1,
        τC = 0,
        condition_on_target = true)

A test for *time-lagged directional influence* between variables
based on the predictive asymmetry principle (Haaga et al., in prep).

## Arguments

- `measure` is some association measure, e.g. [`CMIShannon`](@ref).
- `est` specifies an estimator of `measure`. Defaults to `nothing`. Some measures,
    however, like [`CMIShannon`](@ref), requires an estimator (e.g. [`FPVP`](@ref)).

## Keyword arguments

- **`ηT`**: A set of positive prediction lags for the target variable. Must contain
    strictly positive integers.
- **`τS`**: Embedding delay(s) for the source variable(s). Must contain non-negative integers.
- **`τT`**: Embedding delay(s) for the target variable(s). Must contain non-negative integers.
    This parameter is ignored if `condition_on_target == false`.
- **`τC`**: Embedding delay(s) for the conditional variable(s). Must contain non-negative
    integers. Only applies if three input variables are provided.
- **`condition_on_target`**: If `condition_on_target == true`, then the target variable
    is also embedded (using `τT`) and included in the condition; this requires
    `measure` to be compatible with multivariate inputs.
    If `condition_on_target == false`, then no conditioning is performed and `measure` must
    be compatible with bivariate inputs (e.g. [`MIShannon`](@ref)).

## Input data

Input data must be unlagged univariate timeseries or unlagged multivariate
[`StateSpaceSet`](@ref) timeseries.
"""
Base.@kwdef struct AsymmetryTest{M, E, N, TS, TT, TC, K, B} <: IndependenceTest{M}
    measure::M
    est::E = nothing
    ηT::N = 1:12
    τS::TS = 0
    τT::TT = 0
    τC::TC = 0
    k::K = 1.0
    condition_on_target::B = true

    function AsymmetryTest(measure::M, est::E = nothing;
            ηT::N = 1:12, τS::TS = 0, τT::TT = 0, τC::TC = 0,
            k::K = 1.0, condition_on_target::B = true) where {M, E, N, TS, TT, TC, K, B}
        new{M, E, N, TS, TT, TC, K, B}(measure, est, ηT, τS, τT, τC, k, condition_on_target)
    end
end

"""
    AsymmetryTestResult(Δ, pvalue)

Holds the result of a [`SurrogateTest`](@ref). `Δ` is the asymmetry distribution.
`pvalue` is the right-sided `p`-value for the test.
"""
struct AsymmetryTestResult{D, F, B, P, T} <: IndependenceTestResult
    n_vars::Int # 2 vars = pairwise, 3 vars = conditional
    Δ::D
    fws::F
    bws::B
    pvalue::P
    test::T
end
pvalue(r::AsymmetryTestResult) = r.pvalue

function Base.show(io::IO, test::AsymmetryTestResult)
    print(io,
        """\
        `AsymmetryTest` for independence
        $(null_hypothesis_text(test))
        median(Δ):   $(median(test.Δ))
        Δ:\t   $(test.Δ)
        fws:\t   $(test.fws)
        bws:\t   $(test.bws)
        $(pvalue_text_summary(test))
        """
        )
end

function independence(test::AsymmetryTest, x, y)
    (; measure, est, ηT, τS, τT, τC, k, condition_on_target) = test
    N = length(x)

    X⁻, X⁺ = lag_for_asymmetry(x, τS)

    if condition_on_target
        maxlag = maximum(abs.([ηT...; τT...; τS...]))
        Y⁻, Y⁺ = lag_for_asymmetry(y, τT)
    else
        maxlag = maximum(abs.([ηT...; τS...]))
    end

    # Account for circularly shifting the data.
    idxs = (maxlag + 1):(N - maxlag)

    Δ = zeros(length(ηT))
    fws = zeros(length(ηT))
    bws = zeros(length(ηT))
    for i in ηT
        Yⁿ⁻, Yⁿ⁺ = lag_for_asymmetry(y, abs(i))
        # TODO: make a version that doesn't use conditioning.
        if condition_on_target
            fw = @views dispatch(measure, est, X⁻[idxs], Yⁿ⁺[idxs], Y⁻[idxs])
            bw = @views dispatch(measure, est, X⁺[idxs], Yⁿ⁻[idxs], Y⁺[idxs])
        else
            fw = @views dispatch(measure, est, X⁻[idxs], Yⁿ⁺[idxs])
            bw = @views dispatch(measure, est, X⁺[idxs], Yⁿ⁻[idxs])
        end
        Δ[i] = fw - bw
        fws[i] = fw
        bws[i] = bw
    end
    Δ = scale_exponentially(Δ; k)

    test = SignedRankTest(Δ)
    p = pvalue(test, tail = :right)
    return AsymmetryTestResult(2, Δ, fws, bws, p, test)
end

function independence(test::AsymmetryTest, x, y, z)
    (; measure, est, ηT, τS, τT, τC, k, condition_on_target) = test
    N = length(x)

    X⁻, X⁺ = lag_for_asymmetry(x, τS)
    Z⁻, Z⁺ = lag_for_asymmetry(z, τC)

    if condition_on_target
        Y⁻, Y⁺ = lag_for_asymmetry(y, τT)
        C⁻ = hcat(Y⁻, Z⁻)
        C⁺ = hcat(Y⁺, Z⁺)
        maxlag = maximum(abs.([ηT...; τS...; τT...; τC...]))
    else
        C⁻ = Z⁻
        C⁺ = Z⁺
        maxlag = maximum(abs.([ηT...; τS...; τC...]))
    end

    # Account for circularly shifting the data.
    idxs = (maxlag + 1):(N - maxlag)

    Δ = zeros(length(ηT))
    fws = zeros(length(ηT))
    bws = zeros(length(ηT))

    for i in ηT
        Yⁿ⁻, Yⁿ⁺ = lag_for_asymmetry(y, abs(i))
        fw = @views dispatch(measure, est, X⁻[idxs], Yⁿ⁺[idxs], C⁻[idxs])
        bw = @views dispatch(measure, est, X⁺[idxs], Yⁿ⁻[idxs], C⁺[idxs])
        Δ[i] = fw - bw
        fws[i] = fw
        bws[i] = bw
    end
    Δ = scale_exponentially(Δ; k)
    test = SignedRankTest(Δ)
    p = pvalue(test, tail = :right)
    return AsymmetryTestResult(3, Δ, fws, bws, p, test)
end

function scale_exponentially(Δ; k = 0.5)
    diffs = [1 - i for i in eachindex(Δ)]
    wts = exp.(k .* diffs)
    return Δ .* wts
end

"""
    lag_for_asymmetry(x, τs) → X⁻, X⁺

Given `x`, an unlagged univariate timeseries or unlagged multivariate
[`StateSpaceSet`](@ref) timeseries, return circularly shifted
timeseries `X⁻` and `X⁺`.

If `τs` is an integer, then `x` is circularly shifted by `τs` units, both positively and
negatively. If `τs` is an iterable of many integers, then `X⁻` contains
a matrix where the `i`-th column is `x` circularly shifted by `τs[i]`.

Note: this operation yields `missing` values at the start/end of the
`X⁻` and `X⁺` that may have to be discarded in later operations.
"""
function lag_for_asymmetry(x, τs)
    X⁻ = hcat((ShiftedArrays.circshift(x, i) for i in τs)...)
    X⁺ = hcat((ShiftedArrays.circshift(x, -i) for i in τs)...)
    return Dataset(X⁻), Dataset(X⁺)
end

function lag_for_asymmetry(x::AbstractStateSpaceSet, τs)
    X = Matrix(x)
    X⁻ = hcat((ShiftedArrays.circshift(X, i) for i in τs)...)
    X⁺ = hcat((ShiftedArrays.circshift(X, -i) for i in τs)...)
    return Dataset(X⁻), Dataset(X⁺)
end

# Internal method used time series causal graph inference.
function lag_for_asymmetry(x::AbstractStateSpaceSet, τs, js)
    X = columns(x)
    X⁻ = hcat([ShiftedArrays.circshift(X[j], τ) for (j, τ) in zip(τs, js)]...)
    X⁺ = hcat((ShiftedArrays.circshift(X[j], -τ) for (j, τ) in zip(τs, js))...)
    return Dataset(X⁻), Dataset(X⁺)
end

function lag_for_asymmetry(x, τs, js)
    X⁻ = hcat([ShiftedArrays.circshift(x[j], τ) for (j, τ) in zip(τs, js)]...)
    X⁺ = hcat((ShiftedArrays.circshift(x[j], -τ) for (j, τ) in zip(τs, js))...)
    return Dataset(X⁻), Dataset(X⁺)
end
