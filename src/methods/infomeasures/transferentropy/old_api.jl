
"""
transferentropy([e::EntropyDefinition], s, t, [c,] est;
    τT = -1, τS = -1, η𝒯 = 1, dT = 1, dS = 1, d𝒯 = 1, [τC = -1, dC = 1])

Estimate transfer entropy[^Schreiber2000] from source `s` to target `t`, ``TE^{q}(s \\to t)``, using the
provided entropy/probability estimator `est` with logarithms to the given `base`. Optionally, condition
on `c` and estimate the conditional transfer entropy ``TE^{q}(s \\to t | c)``.

The input timeseries `s`, `t`, and `c` must be equal-length real-valued vectors.
The entropy (first argument) is optional; if not given, `Shannon()` is used instead.

!!! note "Types of transfer entropy"
The first argument `e` controls the type of transfer entropy, and `e` can be any
`Entropy` from Entropies.jl. The default, `Shannon()`, computes the regular Shannon
transfer entropy. Use `Renyi()` to compute generalized Rényi transfer
entropy[^Jizba2012].

All possible estimators that can be used are described in the online documentation.

## Keyword Arguments

Keyword arguments tune the embedding that will be done to each of the timeseries
(with more details following below).
In short, the embedding lags `τT`, `τS`, `τC` must be zero or negative, the
prediction lag `η𝒯` must be positive, and the embedding dimensions `dT`, `dS`, `dC`, `d𝒯`
must be greater than or equal to 1. Thus, the convention is to use negative lags to
indicate embedding delays for past state vectors (for the ``T``, ``S`` and ``C`` marginals,
detailed below), and positive lags to indicate embedding delays for future state vectors
(for the ``\\mathcal T`` marginal, also detailed below).

The default behaviour is to use scalar timeseries for past state vectors
(in that case, the `τT`, `τS` or `τC` does not affect the analysis).

## Description

### Transfer entropy on scalar time series

Transfer entropy[^Schreiber2000] between two simultaneously measured scalar time series ``s(n)`` and ``t(n)``,
``s(n) = \\{ s_1, s_2, \\ldots, s_N \\} `` and ``t(n) = \\{ t_1, t_2, \\ldots, t_N \\} ``, is
is defined as

```math
TE(s \\to t) = \\sum_i p(s_i, t_i, t_{i+\\eta}) \\log \\left( \\dfrac{p(t_{i+\\eta} | t_i, s_i)}{p(t_{i+\\eta} | t_i)} \\right)
```

### Transfer entropy on generalized embeddings

By defining the vector-valued time series, it is possible to include more than one
historical/future value for each marginal (see 'Uniform vs. non-uniform embeddings' below
for embedding details):

- ``\\mathcal{T}^{(d_{\\mathcal T}, \\eta_{\\mathcal T})} = \\{t_i^{(d_{\\mathcal T}, \\eta_{\\mathcal T})} \\}_{i=1}^{N}``,
- ``T^{(d_T, \\tau_T)} = \\{t_i^{(d_T, \\tau_T)} \\}_{i=1}^{N}``,
- ``S^{(d_S, \\tau_S)} = \\{s_i^{(d_T, \\tau_T)} \\}_{i=1}^{N}``,  and
- ``C^{(d_C, \\tau_C)} = \\{s_i^{(d_C, \\tau_C)} \\}_{i=1}^{N}``.

The non-conditioned generalized and conditioned generalized forms of the transfer entropy are then

```math
TE(s \\to t) = \\sum_i p(S,T, \\mathcal{T}) \\log \\left( \\dfrac{p(\\mathcal{T} | T, S)}{p(\\mathcal{T} | T)} \\right)
```

```math
TE(s \\to t | c) = \\sum_i p(S,T, \\mathcal{T}, C) \\log \\left( \\dfrac{p(\\mathcal{T} | T, S, C)}{p(\\mathcal{T} | T, C)} \\right)
```

### Uniform vs. non-uniform embeddings

The `N` state vectors for each marginal are either

- uniform, of the form ``x_{i}^{(d, \\omega)} = (x_i, x_{i+\\omega}, x_{i+2\\omega}, \\ldots x_{i+(d - 1)\\omega})``,
with equally spaced state vector entries. *Note: When constructing marginals for ``T``, ``S`` and ``C``,
we need ``\\omega \\leq 0`` to get present/past values, while ``\\omega > 0`` is necessary to get future states
when constructing ``\\mathcal{T}``.*
- non-uniform, of the form ``x_{i}^{(d, \\omega)} = (x_i, x_{i+\\omega_1}, x_{i+\\omega_2}, \\ldots x_{i+\\omega_{d}})``,
with non-equally spaced state vector entries ``\\omega_1, \\omega_2, \\ldots, \\omega_{d}``,
which can be freely chosen. *Note: When constructing marginals for ``T``, ``S`` and ``C``,
we need ``\\omega_i \\leq 0`` for all ``\\omega_i`` to get present/past values, while ``\\omega_i > 0`` for all ``\\omega_i``
is necessary to get future states when constructing ``\\mathcal{T}``.*

In practice, the `dT`-dimensional, `dS`-dimensional and `dC`-dimensional state vectors
comprising ``T``, ``S`` and ``C`` are constructed with embedding lags `τT`,
`τS`, and `τC`, respectively. The `d𝒯`-dimensional future states ``\\mathcal{T}^{(d_{\\mathcal T}, \\eta_{\\mathcal T})}``
are constructed with prediction lag `η𝒯` (i.e. predictions go from present/past states to
future states spanning a maximum of `d𝒯*η𝒯` time steps).
*Note: in Schreiber's paper, only the historical states are defined as
potentially higher-dimensional, while the future states are always scalar.*

### Estimation

Transfer entropy is here estimated by rewriting the above expressions as a sum of marginal
entropies, and extending the definitions above to use Rényi generalized entropies of order
`q` as

```math
TE^{q}(s \\to t) = H^{q}(\\mathcal T, T) + H^{q}(T, S) - H^{q}(T) - H^{q}(\\mathcal T, T, S),
```

```math
TE^{q}(s \\to t | c) = H^{q}(\\mathcal T, T, C) + H^{q}(T, S, C) - H^{q}(T, C) - H^{q}(\\mathcal T, T, S, C),
```

where ``H^{q}(\\cdot)`` is the generalized Rényi entropy of order ``q``. This is equivalent
to the Rényi transfer entropy implementation in Jizba et al. (2012)[^Jizba2012].

## Examples

Default estimation (scalar marginals):

```julia
# Symbolic estimator, motifs of length 4, uniform delay vectors with lag 1
est = SymbolicPermutation(m = 4, τ = 1)

x, y = rand(100), rand(100)
transferentropy(x, y, est)
```

Increasing the dimensionality of the ``T`` marginal (present/past states of the target
variable):

```julia
# Binning-based estimator
est = VisitationFrequency(RectangularBinning(4))
x, y = rand(100), rand(100)

# Uniform delay vectors when `τT` is an integer (see explanation above)
# Here t_{i}^{(dT, τT)} = (t_i, t_{i+τ}, t_{i+2τ}, \\ldots t_{i+(dT-1)τ})
# = (t_i, t_{i-2}, t_{i-4}, \\ldots t_{i-6τ}), so we need zero/negative values for `τT`.
transferentropy(x, y, est, dT = 4, τT = -2)

# Non-uniform delay vectors when `τT` is a vector of integers
# Here t_{i}^{(dT, τT)} = (t_i, t_{i+τ_{1}}, t_{i+τ_{2}}, \\ldots t_{i+τ_{dT}})
# = (t_i, t_{i-7}, t_{i-25}), so we need zero/negative values for `τT`.
transferentropy(x, y, est, dT = 3, τT = [0, -7, -25])
```

Logarithm bases and the order of the Rényi entropy can also be tuned:

```julia
x, y = rand(100), rand(100)
est = NaiveKernel(0.3)
transferentropy(x, y, est, base = MathConstants.e, q = 2) # TE in nats, order-2 Rényi entropy
```

[^Schreiber2000]: Schreiber, T. (2000). Measuring information transfer. Physical review letters, 85(2), 461.
[^Jizba2012]: Jizba, P., Kleinert, H., & Shefaat, M. (2012). Rényi’s information transfer between financial time series. Physica A: Statistical Mechanics and its Applications, 391(10), 2971-2989.
"""
# function transferentropy end
# function transferentropy! end


function transferentropy(e::EntropyDefinition, est::ProbabilitiesEstimator,
    args...; base = 2, params = params = TransferEntropy(), kwargs...)

    transferentropy(Shannon(; base), est, args...; kwargs...)
end
function transferentropy(est::TransferDifferentialEntropyEstimator, args...; base = 2, kwargs...)
    transferentropy(Shannon(; base), est, args...; kwargs...)
end



"""
    transferentropy([e::EntropyDefinition,] est::ProbabilitiesEstimator, s, t, [c];
        τT = -1, τS = -1, η𝒯 = 1, dT = 1, dS = 1, d𝒯 = 1, [τC = -1, dC = 1])
    transferentropy(]e::EntropyDefinition,] est::TransferDifferentialEntropyEstimator, s, t, [c]; kwargs...)

Estimate transfer entropy from `s` to `t`, optionally conditioning on `c`,
by a sum of marginal entropies.

# Arguments

- `e::EntropyDefinition`. The type of entropy to compute. Optional. Defaults to `Shannon(; base = 2)`
    if not specified. The unit of the transfer entropy is controlled by the logarithm base,
    which must be provided as a keyword to `e`.
- `est`. Either a [`ProbabilitiesEstimator`](@ref) or a [`TransferDifferentialEntropyEstimator`](@ref),
    which controls how probabilities and entropies are estimated for each marginal.
- `s::AbstractVector`: The source time series
- `t::AbstractVector`: The target time series
- `c::AbstractVector`: An optional time series to condition on.

See also: [`Entropy`](@ref), [`ProbabilitiesEstimator`](@ref).
"""
function transferentropy(e::EntropyDefinition, est::ProbabilitiesEstimator,
        args...; params = TransferEntropy(), kwargs...)
    joint, ST, T𝒯, T = get_marginals(params, args...; emb = EmbeddingTE(; kwargs...))
    from_marginals(TE(), e, est, joint, ST, T𝒯, T)
end

"""
    transferentropy(e::DifferentialEntropyEstimator, s, t, [c];
        τT = -1, τS = -1, η𝒯 = 1, dT = 1, dS = 1, d𝒯 = 1, [τC = -1, dC = 1])

Estimate Shannon transfer entropy from `s` to `t`, optionally conditioning on `c`,
by a sum of marginal entropies, using an [`DifferentialEntropyEstimator`](@ref) estimator from
[Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl).

These methods estimate entropies using some procedure that doesn't explicitly construct a
probability distribution.
"""
function transferentropy(e::DifferentialEntropyEstimator, args...; kwargs...)
    joint, ST, T𝒯, T = get_marginals(TE(), args...; emb = EmbeddingTE(; kwargs...))
    from_marginals(TE(), e, joint, ST, T𝒯, T)
end
