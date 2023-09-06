using StateSpaceSets: StateSpaceSet
using DelayEmbeddings: genembed
import DelayEmbeddings: embed
export PA
export asymmetry

"""
    PA <: CausalityTools.AssociationMeasure
    PA(ηT = 1:5, τS = 1, τC = 1)

The modified predictive asymmetry measure (Haaga et al., in revision).

!!! note
    This is an experimental measure. It is part of an ongoing paper submission revision,
    but is provided here for convenience.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for pairwise
    or conditional directional dependence.
- Use with [`asymmetry`](@ref) to compute the raw asymmetry distribution.

## Keyword arguments

- `ηT`. The prediction lags for the target variable.
- `τS`. The embedding delay(s) for the source variable.
- `τC`. The embedding delay(s) for the conditional variable(s).

All parameters are given as a single integer or multiple integers greater than zero.

## Compatible estimators

`PA`/[`asymmetry`](@ref) uses [`condmutualinfo`](@ref) under the hood. Any estimator that
can be used for [`ConditionalMutualInformation`](@ref) can therefore, in principle,
be used with the predictive asymmetry. We recommend to use [`FPVP`](@ref), or one
of the other dedicated conditional mutual information estimators.

| Estimator                        | Type                                            | Principle           | Pairwise | Conditional |
| -------------------------------- | ----------------------------------------------- | ------------------- | :------: | :---------: |
| [`UniqueElements`](@ref)         | [`ProbabilitiesEstimator`](@ref)                | Frequencies         |    ✓    |     ✓      |
| [`ValueHistogram`](@ref)         | [`ProbabilitiesEstimator`](@ref)                | Binning (histogram) |    ✓    |     ✓      |
| [`Dispersion`](@ref)             | [`ProbabilitiesEstimator`](@ref)                | Dispersion patterns |    ✓    |     ✓      |
| [`Kraskov`](@ref)                | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |    ✓    |     ✓      |
| [`Zhu`](@ref)                    | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |    ✓    |     ✓      |
| [`ZhuSingh`](@ref)               | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |    ✓    |     ✓      |
| [`Gao`](@ref)                    | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |    ✓    |     ✓      |
| [`Goria`](@ref)                  | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |    ✓    |     ✓      |
| [`Lord`](@ref)                   | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |    ✓    |     ✓      |
| [`LeonenkoProzantoSavani`](@ref) | [`DifferentialEntropyEstimator`](@ref)          | Nearest neighbors   |    ✓    |     ✓      |
| [`GaussanMI`](@ref)              | [`MutualInformationEstimator`](@ref)            | Parametric          |    ✓    |     ✓      |
| [`KSG1`](@ref)                   | [`MutualInformationEstimator`](@ref)            | Continuous          |    ✓    |     ✓      |
| [`KSG2`](@ref)                   | [`MutualInformationEstimator`](@ref)            | Continuous          |    ✓    |     ✓      |
| [`GaoKannanOhViswanath`](@ref)   | [`MutualInformationEstimator`](@ref)            | Mixed               |    ✓    |     ✓      |
| [`GaoOhViswanath`](@ref)         | [`MutualInformationEstimator`](@ref)            | Continuous          |    ✓    |     ✓      |
| [`FPVP`](@ref)                   | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |    ✓    |     ✓      |
| [`MesnerShalizi`](@ref)          | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |    ✓    |     ✓      |
| [`Rahimzamani`](@ref)            | [`ConditionalMutualInformationEstimator`](@ref) | Nearest neighbors   |    ✓    |     ✓      |

## Examples

- [Computing the asymmetry distribution](@ref examples_pa_asymmetry_dist).
"""
struct PA{N, TS, TC} <: DirectedAssociationMeasure
    ηT::N # multiple positive, unique integers (or some instance of a delay optimiser)
    τS::TS # integer or range of positive integers (or some instance of a delay optimiser)
    τC::TC# integer or range of positive integers (or some instance of a delay optimiser)

    function PA(; ηT::N = 1:5, τS::TS = 1, τC::TC = 1) where {N, TS, TC}
       new{N, TS, TC}(ηT, τS, τC)
    end
end

# Optimize embedding lags. The following method just returns the given lags,
# whatever form they may be in. Other methods of `optimize_embedding`
# run some optimization procedure.
optimize_embedding(τs, x) = τs

include("embeddings.jl")

"""
    asymmetry(measure::PA, est, x, y, [z]) → ΔA

Compute the predictive asymmetry ([`PA`](@ref)) from `x` to `y`, conditioned on `z` if
given. Returns the distribution of asymmetry values.

## Compatible estimators

The docstring for [`PA`](@ref) lists compatible estimators.
"""
function asymmetry(args...; kwargs...)
    return estimate(args...; kwargs...)
end

const PA_ESTS = Union{
    ProbabilitiesEstimator,
    DifferentialEntropyEstimator,
    MutualInformationEstimator,
    ConditionalMutualInformationEstimator
}

function asymmetry(est::PA_ESTS, args...)
    throw(ArgumentError("A valid measure must be provided as the second argument; do `PA()`."))
end

function estimate(measure::PA, x::AbstractVector...)
    throw(ArgumentError("A valid estimator must be provided as the second argument, try for example `FPVP()`."))
end

function redefine_lags(measure::PA, s)
    if measure.τS isa Pecuzal
        τS = optimize_embedding(measure.τS, s)
    else
        τS = measure.τS
    end
    return PA(; ηT = measure.ηT, τS = τS, τC = measure.τC)
end

function redefine_lags(measure::PA, s, c)
    if measure.τS isa Pecuzal
        τS = optimize_embedding(measure.τS, s)
    else
        τS = measure.τS
    end

    if measure.τC isa Pecuzal
        τC = optimize_embedding(measure.τC, c)
    else
        τC = measure.τC
    end
    return PA(; ηT = measure.ηT, τS = τS, τC = τC)
end

function embed_for_pa(measure::PA, x, y)
    (; ηT, τS, τC) = measure

    ηT = sort(unique(ηT))
    τS = sort(unique([τS; ηT]))
    all(ηT .> 0) || throw(ArgumentError("Lags must be positive integers. Got ηT=$ηT."))
    all(τS .> 0) || throw(ArgumentError("Lags must be positive integers. Got τS=$τS."))
    all(τC .> 0) || throw(ArgumentError("Lags must be positive integers. Got τC=$τC."))

    nη = length(ηT)
    nτx = length(τS)
    jη = repeat([2], nη)
    js_x = repeat([1], nτx)
    js = [jη...; 2; jη...;     js_x...; 1; js_x...]
    τs = [ηT...; 0; .-(ηT)...; τS...; 0; .-(τS)...]
    joint = genembed(StateSpaceSet(x, y), τs, js)
    T⁺ = joint[:, 1:nη]
    T⁰ = joint[:, nη+1:nη+1]
    T⁻ = joint[:, nη+2:2nη+1]
    iₛ = 2nη+2
    S⁺ = joint[:, iₛ:(iₛ+nτx-1)]
    S⁰ = joint[:, (iₛ+nτx):(iₛ+nτx)]
    S⁻ = joint[:, (iₛ+nτx+1):size(joint, 2)]

    return T⁺, T⁰, T⁻, S⁺, S⁰, S⁻
end

function embed_for_pa(measure::PA, x, y, z)
    (; ηT, τS, τC) = measure
    ηT = sort(unique(ηT))
    τS = sort(unique([τS; ηT]))
    τC = sort(unique([τC; ηT]))
    all(ηT .> 0) || throw(ArgumentError("Lags must be positive integers. Got ηT=$ηT."))
    all(τS .> 0) || throw(ArgumentError("Lags must be positive integers. Got τS=$τS."))
    all(τC .> 0) || throw(ArgumentError("Lags must be positive integers. Got τC=$τC."))
    nη = length(ηT)
    nτx = length(τS)
    nτz = length(τC)

    jη = repeat([2], nη)
    js_x = repeat([1], nτx)
    js_z = repeat([3], nτz)

    js = [jη...; 2; jη...;     js_x...; 1; js_x...;   js_z...; 1; js_z...]
    τs = [ηT...; 0; .-(ηT)...; τS...;   0; .-(τS)...; τC...;   0; .-(τC)...]
    joint = genembed(StateSpaceSet(x, y, z), τs, js)
    T⁺ = joint[:, 1:nη]
    T⁰ = joint[:, nη+1:nη+1]
    T⁻ = joint[:, nη+2:2nη+1]
    iₛ = 2nη+2
    S⁺ = joint[:, iₛ:(iₛ+nτx-1)]
    S⁰ = joint[:, (iₛ+nτx):(iₛ+nτx)]
    S⁻ = joint[:, (iₛ+nτx+1):(iₛ+2nτx)]
    ic = iₛ+2nτx+1
    C⁺ = joint[:, ic:(ic+nτz-1)]
    C⁰ = joint[:, (ic+nτz):(ic+nτz)]
    C⁻ = joint[:, (ic+nτz+1):(ic+2nτz)]

    return T⁺, T⁰, T⁻, S⁺, S⁰, S⁻, C⁺, C⁰, C⁻
end

as_vector(x::Union{AbstractVector{T}}) where T = vec(x)
as_vector(x::AbstractStateSpaceSet{1, T}) where T = x[:, 1]

function estimate(measure::PA, est::PA_ESTS, x, y)
    X = as_vector(x)
    Y = as_vector(y)
    measure = redefine_lags(measure, x)
    (; ηT, τS, τC) = measure
    T⁺, T⁰, T⁻, S⁺, S⁰, S⁻ = embed_for_pa(measure, X, Y)

    N = length(T⁺)
    τŜ = sort(unique([τS; ηT]))
    S⁺⁰ = StateSpaceSet(S⁺[:, findall(x -> x ∈ τS, τŜ)], S⁰)
    S⁻⁰ = StateSpaceSet(S⁻[:, findall(x -> x ∈ τS, τŜ)], S⁰)
    ΔA = zeros(length(ηT))
    fw = zeros(length(ηT))
    bw = zeros(length(ηT))
    T⁺⁰ = MVector{2}.(StateSpaceSet(T⁰, T⁰).data)
    T⁻⁰ = MVector{2}.(StateSpaceSet(T⁰, T⁰).data)
    Tⁿ⁺ = MVector{1}.(StateSpaceSet(T⁰).data)
    Tⁿ⁻ = MVector{1}.(StateSpaceSet(T⁰).data)
    for (i, η) in enumerate(sort(ηT))
        fill_target!(Tⁿ⁺, T⁺, i)
        fill_target!(Tⁿ⁻, T⁻, i)
        fill_target!(T⁺⁰, T⁺, i)
        fill_target!(T⁻⁰, T⁻, i)
        fw[i] = condmutualinfo(CMIShannon(), est, Tⁿ⁺, S⁻⁰, T⁻⁰)
        bw[i] = condmutualinfo(CMIShannon(), est, Tⁿ⁻, S⁺⁰, T⁺⁰)
        ΔA[i] = fw[i] - bw[i]
    end

    return ΔA
end

function estimate(measure::PA, est::PA_ESTIMATORS, x, y, z)
    measure = redefine_lags(measure, x, z)
    (; ηT, τS, τC) = measure
    X = as_vector(x)
    Y = as_vector(y)
    Z = as_vector(z)
    T⁺, T⁰, T⁻, S⁺, S⁰, S⁻, C⁺, C⁰, C⁻ = embed_for_pa(measure, X, Y, Z)

    N = length(T⁺)
    τŜ = sort(unique([τS; ηT]))
    S⁺⁰ = StateSpaceSet(S⁺[:, findall(x -> x ∈ τS, τŜ)], S⁰)
    S⁻⁰ = StateSpaceSet(S⁻[:, findall(x -> x ∈ τS, τŜ)], S⁰)
    τĈ = sort(unique([τC; ηT]))
    C⁺⁰ = StateSpaceSet(C⁺[:, findall(x -> x ∈ τC, τĈ)], C⁰)
    C⁻⁰ = StateSpaceSet(C⁻[:, findall(x -> x ∈ τC, τĈ)], C⁰)
    ΔA = zeros(length(ηT))
    fw = zeros(length(ηT))
    bw = zeros(length(ηT))
    Tⁿ⁺⁰C⁺⁰ = MVector{2 + dimension(C⁻⁰)}.(StateSpaceSet(T⁰, T⁰, C⁻⁰).data)
    Tⁿ⁻⁰C⁻⁰ = MVector{2 + dimension(C⁺⁰)}.(StateSpaceSet(T⁰, T⁰, C⁺⁰).data)
    Tⁿ⁺ = MVector{1}.(T⁰.data)
    Tⁿ⁻ = MVector{1}.(T⁰.data)

    for (i, η) in enumerate(sort(ηT))
        fill_target!(Tⁿ⁺, T⁺, i)
        fill_target!(Tⁿ⁻, T⁻, i)
        fill_target!(Tⁿ⁺⁰C⁺⁰, T⁺, i)
        fill_target!(Tⁿ⁻⁰C⁻⁰, T⁺, i)
        fw[i] = condmutualinfo(CMIShannon(), est, Tⁿ⁺, S⁻⁰, Tⁿ⁻⁰C⁻⁰)
        bw[i] = condmutualinfo(CMIShannon(), est, Tⁿ⁻, S⁺⁰, Tⁿ⁺⁰C⁺⁰)
        ΔA[i] = fw[i] - bw[i]
    end

    return ΔA
end

# The same also applies for T⁻⁰ and T⁻, just exchange variables
function fill_target!(T⁺⁰, T⁺, i, k::Int = 1)
    for j in eachindex(T⁺⁰)
        T⁺⁰[j][k] = T⁺[j, i]
    end
end
