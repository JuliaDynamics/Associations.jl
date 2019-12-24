import TransferEntropy: transferentropy

"""
    transferentropy(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int,
        binning_scheme::Union{RectangularBinning, Vector{RectangularBinning}}; 
        summary_statistic::Function = StatsBase.mean,
        η::Int = 1, τ::Int = 1, 
        estimator::BinningTransferEntropyEstimator = VisitationFrequency(b = 2)) -> FLoat64

Convenience function for computing transfer entropy from a `source` time series to a 
`response` time series. 

Transfer entropy is computed using the provided `estimator` over 
a `k + l + m`-dimensional [generalised delay reconstruction](@ref custom_delay_reconstruction) 
of the input data. The embedding delay used for the reconstructions is `τ`, 
and the prediction lag is `η`. Transfer entropy is computed using the following conventions 
for the [generalised delay reconstructions](@ref custom_delay_reconstruction).

```math
\\begin{align}
T_f^{(k)} &= \\{(T(t+\\eta_k), \\ldots, T(t+\\eta_2), T(t+\\eta_1)) \\} \\\\
T_{pp}^{(l)} &= \\{ (T(t), T(t-\\tau_1), T(t-\\tau_2), \\ldots, T(t - \\tau_{l - 1})) \\} \\\\
S_{pp}^{(m)} &= \\{ (S(t), S(t-\\tau_1), S(t-\\tau_2), \\ldots, S(t-\\tau_{m - 1})) \\}
\\end{align}
```

where ``T_f`` denotes the `k`-dimensional set of vectors furnishing the future states of ``T``,
``T_{pp}`` denotes the `l`-dimensional set of vectors furnishing the past and present states of ``T``, 
and ``S_{pp}`` denotes the `m`-dimensional set of vectors furnishing the past and present of ``S``. 
``\\eta`` is the prediction lag. This convenience function uses ``\\tau_1`` = `τ`, 
``\\tau_2`` = `2*τ`, ``\\tau_3`` = `3*τ`, and so on.

Combined, we get the generalised embedding ``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)})``.
TE is then computed as 

```math
\\begin{align}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}) \\log_{b}{\\left(\\frac{P(T_f | T_{pp}, S_{pp})}{P(T_f | T_{pp})}\\right)}
\\end{align}
```

The probabilities needed for computing ``TE_{S \\rightarrow T}`` are estimated 
over a discretization of the state space reconstruction, and this discretization is 
dictated by the provided `binning_scheme`(s). If there is more than one binning scheme,
then the transfer entropy is summarised using `summary_statistic` over the partitions, 
and a single value for the transfer entropy is returned.

## Arguments

- **`source`**: The source data series.
- **`target`**: The target data series.
- **`k`**: The dimension of the ``T_{f}`` component of the embedding. 
- **`l`**: The dimension of the ``T_{pp}`` component of the embedding. 
- **`m`**: The dimension of the ``S_{pp}`` component of the embedding. 
- **`binning_scheme`**: The binning scheme(s) used to construct the partitions
    over which TE is computed. Must be either one or several instances of 
    `RectangularBinning`s (provided as a vector). TE is computed for each 
    of the resulting partitions.

## Keyword arguments

- **`τ::Int = 1`**: The embedding lag (lags the ``T_{pp}`` component of the embedding). 
- **`η::Int = 1`**: The prediction lag. 
- **`estimator::BinningTransferEntropyEstimator = VisitationFrequency()`**: The 
    transfer entropy estimator to use.

## Returns 

A single value for the transfer entropy.
"""
function transferentropy(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int,
        binning_scheme::Union{RectangularBinning, Vector{RectangularBinning}}; 
        summary_statistic::Function = StatsBase.mean,
        η = 1, τ = 1, 
        estimator = VisitationFrequency(b = 2))

    k + l + m >= 3 || throw(ArgumentError("`dim = k + l + m` must be 3 or higher for regular TE"))

    pts, vars = te_embed(source, response, k, l, m, η = η, τ = τ)
    
    # Compute TE over different partitions
    # ====================================
    bs = binning_scheme isa Vector{RectangularBinning} ? binning_scheme : [binning_scheme]
    tes = map(binscheme -> transferentropy(pts, vars, binscheme, estimator), bs)

    return summary_statistic(tes)
end

"""
    transferentropy(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        cond::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int,
        binning_scheme::Union{RectangularBinning, Vector{RectangularBinning}}; 
        summary_statistic::Function = StatsBase.mean,
        η::Int = 1, τ::Int = 1, 
        estimator::BinningTransferEntropyEstimator = VisitationFrequency(b = 2)) -> FLoat64

Convenience function for computing transfer entropy from a `source` time series to a 
`response` time series, conditioned on a third time series `cond`.

Transfer entropy is computed using the provided `estimator` over 
a `k + l + m + n`-dimensional [generalised delay reconstruction](@ref custom_delay_reconstruction) 
of the input data. The embedding delay used for the reconstructions is `τ`, 
and the prediction lag is `η`. Transfer entropy is computed using the following conventions 
for the [generalised delay reconstructions](@ref custom_delay_reconstruction).

```math
\\begin{align}
T_f^{(k)} &= \\{(T(t+\\eta_k), \\ldots, T(t+\\eta_2), T(t+\\eta_1)) \\} \\\\
T_{pp}^{(l)} &= \\{ (T(t), T(t-\\tau_1), T(t-\\tau_2), \\ldots, T(t - \\tau_{l - 1})) \\} \\\\
S_{pp}^{(m)} &= \\{ (S(t), S(t-\\tau_1), S(t-\\tau_2), \\ldots, S(t-\\tau_{m - 1})) \\} \\\\
C_{pp}^{(n)} &= \\{ (C(t), C(t-\\tau_1), C(t-\\tau_2), \\ldots, C(t-\\tau_{n - 1})) \\}
\\end{align}
```

where ``T_f`` denotes the `k`-dimensional set of vectors furnishing the future states of ``T``,
``T_{pp}`` denotes the `l`-dimensional set of vectors furnishing the past and present states of ``T``, 
``S_{pp}`` denotes the `m`-dimensional set of vectors furnishing the past and present of ``S``, 
and ``C_{pp}`` denotes the `n`-dimensional set of vectors furnishing the past and present of ``C``.
``\\eta`` is the prediction lag. This convenience function uses ``\\tau_1`` = `τ`, 
``\\tau_2`` = `2*τ`, ``\\tau_3`` = `3*τ`, and so on. 

Combined, we get the generalised embedding  ``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)}, C_{pp}^{(n)})``.
TE is then computed as 

```math
\\begin{align}
TE_{S \\rightarrow T | C} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}, C_{pp}) \\log_{b}{\\left(\\frac{P(T_f | T_{pp}, S_{pp}, C_{pp})}{P(T_f | T_{pp}, C_{pp})}\\right)}
\\end{align}
```

The probabilities needed for computing ``TE_{S \\rightarrow T | C}`` are estimated 
over a discretization of the state space reconstruction, and this discretization is 
dictated by the provided `binning_scheme`(s). If there is more than one binning scheme,
then the transfer entropy is summarised using `summary_statistic` over the partitions, 
and a single value for the transfer entropy is returned.

## Arguments

- **`source`**: The source data series.
- **`target`**: The target data series.
- **`k::Int`**: The dimension of the ``T_{f}`` component of the embedding. 
- **`l::Int`**: The dimension of the ``T_{pp}`` component of the embedding. 
- **`m::Int`**: The dimension of the ``S_{pp}`` component of the embedding. 
- **`n::Int`**: The dimension of the ``C_{pp}`` component of the embedding. 
- **`binning_scheme`**: The binning scheme(s) used to construct the partitions
    over which TE is computed. Must be either one or several instances of 
    `RectangularBinning`s (provided as a vector). TE is computed for each 
    of the resulting partitions.

## Keyword arguments

- **`τ::Int = 1`**: The embedding lag (lags the ``T_{pp}`` component of the embedding). 
- **`η::Int = 1`**: The prediction lag. 
- **`estimator::BinningTransferEntropyEstimator = VisitationFrequency()`**: The 
    transfer entropy estimator to use.

## Returns 

A single value for the transfer entropy.
"""
function transferentropy(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1},
        cond::AbstractArray{<:Real, 1},
        k::Int, l::Int, m::Int, n::Int,
        binning_scheme::Union{RectangularBinning, Vector{RectangularBinning}}; 
        summary_statistic::Function = StatsBase.mean,
        η = 1, τ = 1, 
        estimator = VisitationFrequency(b = 2))

    k + l + m + n >= 4 || throw(ArgumentError("`dim = k + l + m + n` must be 4 or higher for conditional TE"))

    pts, vars = te_embed(source, response, cond, k, l, m, n, τ = τ, η = η)

    # Compute TE over different partitions
    # ====================================
    bs = binning_scheme isa Vector{RectangularBinning} ? binning_scheme : [binning_scheme]
    tes = map(binscheme -> transferentropy(pts, vars, binscheme, estimator), bs)

    return summary_statistic(tes)
end

export transferentropy