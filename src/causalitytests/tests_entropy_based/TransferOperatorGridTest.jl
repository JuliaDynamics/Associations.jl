
import StatsBase

"""
    TransferOperatorGridTest(k::Int = 1, l::Int = 1, m::Int = 1, n::Int = 1, 
        τ::Int = 1, b = 2, estimator::TransferOperatorGrid = TransferOperatorGrid(), 
        binning_summary_statistic::Function = StatsBase.mean,
        binning::RectangularBinning, ηs)

The parameters for a transfer entropy test using the `TransferOperatorGrid` estimator [1].

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

- **`estimator::VisitationFrequency`**: A `VisitationFrequency` estimator instance.

- **`binning_summary_statistic::Function`**: A summary statistic to summarise the 
    transfer entropy values if multiple binnings are provided.


## Estimation of the invariant measure

With a `TransferOperatorGridTest`, the first step is to compute an approximation to the 
[transfer operator](@ref transfer_operator_rectangular)
over the partition elements of a rectangular [discretization](@ref discretization) 
of an appropriate [delay reconstruction](@ref custom_delay_reconstruction) of the time 
series to be analysed. Transfer entropy is then computed from the [invariant 
distribution](@ref invariant_measure_rectangular) arising from the transfer operator. 

## About the delay reconstruction for transfer entropy analysis

Denote the time series for the source process ``S`` as ``S(t)``, and the time series for 
the target process ``T`` as ``T(t)``, and ``C_i(t)`` as the time series for any conditional 
processes ``C_i`` that also may influence ``T``. To compute (conditional) TE, we need a 
generalised embedding [3, 4] incorporating all of these processes.

For convenience, define the state vectors

```math
\\begin{align}
T_f^{(k)} &= \\{(T(t+\\eta_k), \\ldots, T(t+\\eta_2), T(t+\\eta_1))\\}, \\label{eq:Tf} \\\\
T_{pp}^{(l)} &= \\{ (T(t), T(t-\\tau_1), T(t-\\tau_2), \\ldots, T(t - \\tau_{l - 1})) \\\\
S_{pp}^{(m)} &= \\{(S(t), S(t-\\tau_1), S(t-\\tau_2), \\ldots, S(t-\\tau_{m - 1}))\\},\\\\
C_{pp}^{(n)} &= \\{ (C_1(t), C_1(t-\\tau_1), \\ldots,  C_2(t), C_2(t-\\tau_1) \\},
\\end{align}
```

where the state vectors ``T_f^{(k)}`` contain ``k`` future values of the target 
variable, ``T_{pp}^{(l)}`` contain ``l`` present and past values of the target 
variable, ``S_{pp}^{(m)}`` contain ``m`` present and past values of the source 
variable, ``C_{pp}^{(n)}`` contain a total of ``n`` present and past values of any 
conditional variable(s).  Combining all variables, we have the generalised embedding 

```math 
\\begin{align}
\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)}, C_{pp}^{(n)})
\\end{align}
```

with a total embedding dimension of ``k + l + m + n``. 
Here, ``\\tau`` indicates the [embedding](@ref custom_delay_reconstruction) lag 
(in the `VisitationFrequencyTest`, we set ``\\tau_1 = \\tau_2 = \\tau_3 = \\ldots``), and
``\\eta`` indicates the prediction lag (the lag of the influence the source 
has on the target). 


Hence, in the generalised embedding, only ``T_f`` depends on the prediction lag ``\\eta``, 
which is to be determined by the analyst. For transfer entropy analysis, ``\\eta`` is chosen 
to be some positive integer, while for 
[predictive asymmetry analysis](@ref predictive_asymmetry), 
symmetric negative and positive ``\\eta``s are used for computing ``\\mathbb{A}``.  

In terms of this generalised embedding, transfer entropy from a source variable ``S`` to a 
target variable ``T`` with conditioning on variable(s) ``C`` is thus defined as 

```math
\\begin{align}
TE_{S \\rightarrow T|C} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}, C_{pp}) \\log_{2}{\\left( \\frac{P(T_f | T_{pp}, S_{pp}, C_{pp})}{P(T_f | T_{pp}, C_{pp})}\\right)}
\\end{align}
```

Without conditioning, we have 

```math
\\begin{align}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}) \\log_{2}{\\left(\\frac{P(T_f | T_{pp}, S_{pp})}{P(T_f | T_{pp})}\\right)}
\\end{align}
```


## Low-level estimator 

This test uses the `TransferOperatorGrid` estimator on the following low-level method 
under the hood. 

- [`transferentropy(::Any, ::TEVars, ::RectangularBinning, ::TransferEntropyEstimator)`](@ref)

In this estimator, the mapping between variables of the 
[generalised embedding](@ref custom_delay_reconstruction) and the marginals during 
transfer entropy computation is controlled using a [`TEVars`](@ref) 
instance. It is *highly* recommended that you check the documentation for this 
method, because it describes the transfer entropy estimation procedure in detail.

## Notes:

- Use `causality(source, target, params::TransferOperatorGridTest)` for regular 
    transfer entropy analysis. This method uses only the `k`, `l`, `m` and ignores `n` 
    when constructing the delay reconstruction. 

- Use `causality(source, target, cond, params::TransferOperatorGridTest)` for conditional 
    transfer entropy analysis. This method uses the `k`, `l`, `m` *and* `n` when constructing 
    the delay reconstruction.

## Example

```julia
# Prediction lags
ηs = 1:10
binning = RectangularBinning(10)

# Use defaults, binning and prediction lags are required. 
# Note that `binning` and `ηs` are *mandatory* keyword arguments.
TransferOperatorGridTest(binning = binning, ηs = ηs)

# The other keywords can also be adjusted
TransferOperatorGridTest(k = 1, l = 2, binning = binning, ηs = ηs)
```

## References 

1. Diego, David, Kristian Agasøster Haaga, and Bjarte Hannisdal. "Transfer entropy computation 
    using the Perron-Frobenius operator." Physical Review E 99.4 (2019): 042212.
    [https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.042212](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.042212)
"""
Base.@kwdef struct TransferOperatorGridTest <: TransferEntropyCausalityTest
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

    """ The transfer entropy estimator. """
    estimator::TransferOperatorGrid = TransferOperatorGrid()

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
end


function causality(source::AbstractVector{T}, target::AbstractVector{T}, p::TransferOperatorGridTest)  where {T <: Real}
    [p.binning_summary_statistic(
        transferentropy(source, target, p.binning, 
            p.k, p.l, p.m, η = η, τ = p.τ, 
            estimator = p.estimator)) for η in p.ηs]
end

export TransferOperatorGridTest
