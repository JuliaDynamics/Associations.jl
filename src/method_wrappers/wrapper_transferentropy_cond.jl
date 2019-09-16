
"""
    te_cond(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1},
        cond::AbstractArray{<:Real, 1},
        k::Int, l::Int, m::Int, n::Int; 
        η = 1, τ = 1, 
        estimator = VisitationFrequency(), 
        n_subdivs = 1,
        b = 2)

# Conditional TE with default discretization scheme(s)

Calculate transfer entropy from `source` to `response` conditioned on
`cond` using the provided 
`estimator` on a rectangular discretization of a `k + l + m + n`-dimensional 
delay embedding of the input data, using an embedding delay of `τ` across 
all embedding components. `η` is the prediction lag. 

## Arguments
- **`source`**: The source data series.
- **`target`**: The target data series.
- **`cond`**: The conditional data series.
- **`k`**: The dimension of the ``T_{f}`` component of the embedding. 
- **`l`**: The dimension of the ``T_{pp}`` component of the embedding. 
- **`m`**: The dimension of the ``S_{pp}`` component of the embedding. 
- **`n`**: The dimension of the ``C_{pp}`` component of the embedding. 

## Keyword arguments

- **`τ`**: The embedding lag. Default is `τ = 1`.
- **`η`**: The prediction lag. Default is `η = 1`.
- **`estimator`**: The transfer entropy estimator to use. The default 
    is `VisitationFrequency()`.
- **`n_subdivs`**: The number of different partitions of varying coarseness
    to compute TE over, as described below. Default is `n_subdivs = 2`.
    (this way, TE is computed over two separate partitions). Unless `n_subdivs = 0`,
    make sure that `n_subdivs` is the same across analyses if they 
    are to be compared. T TE 
- **`b`**: Base of the logarithm. The default (`b = 2`) gives the TE in bits.

## More about the embedding

To compute transfer entropy, we need an appropriate delay embedding 
of `source` (``S``), `target` (``T``) and `cond` (``C``). For convenience, define 


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
,``S_{pp}`` denotes the `m`-dimensional set of vectors furnishing the past and present of ``S``, 
and ``C_{pp}`` denotes the `n`-dimensional set of vectors furnishing the past and present of ``C``.
``\\eta`` is the prediction lag. This convenience function uses ``\\tau_1`` = `τ`, 
``\\tau_2`` = `2*τ`, ``\\tau_3`` = `3*τ`, and so on.


Combined, we get the generalised embedding ``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)}, C_{pp}^{(n)})``, 
which is discretized as described below.

## More about discretization

To compute TE, we coarse-grain the `k+l+m+n`-dimensional generalised embedding 
``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)}, C_{pp}^{(n)})`` into hyperrectangular boxes.
The magnitude of the TE may be biased by the particular choice of binning scheme,
so we compute TE across a number of different box sizes, determined as follows.

Let ``L`` be the number of observations in ``S`` (and ``T`` and ``C``). The coarsest box size is 
determined by subdiving the i-th coordinate axis into 
``N = ceiling(L^\\frac{1}{k + l + m + n + 1})`` intervals of equal lengths, 
resulting in a box size of ``|max(dim_{i}) - min(dim_{i})|/N``. The next box size 
is given by ``|max(dim_{i}) - min(dim_{i})|/(N+1)``, then 
``|max(dim_{i}) - min(dim_{i})|/(N+2)``, and so on, until the finest box 
size which is given by ``|max(dim_{i}) - min(dim_{i})|/(N+N_{subdivs})``. 

## Transfer entropy computation

TE is then computed as 

```math
\\begin{align}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}, C_{pp}) \\log_{b}{\\left(\\frac{P(T_f | T_{pp}, S_{pp}, C_{pp})}{P(T_f | T_{pp}, C_{pp})}\\right)}
\\end{align}
```
using the provided `estimator` (default = `VisitationFrequency`) for 
each of the discretizations. A vector of the TE estimates for each discretization 
is returned.
"""
function te_cond(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        cond::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int, n::Int; 
        η = 1, τ = 1, 
        estimator::TransferEntropyEstimator = VisitationFrequency(), 
        n_subdivs = 2,
        b = 2)

    k + l + m + n >= 4 || throw(ArgumentError("`dim = k + l + m + n` must be 4 or higher for conditional TE"))

    pts, vars = te_embed(source, response, cond, k, l, m, n, η = η, τ = τ)

    # Determine appropriate binnings from time series length (roughly according to 
    # Krakovska et al. (2018)'s recommendations)
    L = length(source) 
    
    # one subdivision coarser than Krakovska as the coarest partition
    nbins_coarsest = ceil(Int, L^(1/(k + l + m + n + 1))) - 1 
    
    # with the default n_subdivs = 2, we cover both the Krakovska recommendation,
    # as well as one finer bin size and one coarser bin size.
    nbins_finest = nbins_coarsest + n_subdivs 
    binning_scheme = map(n-> RectangularBinning(n), nbins_coarsest:nbins_finest)
    
    # Compute TE over different partitions
    # ====================================
    tes = map(binscheme -> transferentropy(pts, vars, binscheme, estimator, b = b), binning_scheme)

    return tes
end

"""
    te_cond(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1},
        cond::AbstractArray{<:Real, 1},
        k::Int, l::Int, m::Int, n::Int,
        binning_scheme::Vector{RectangularBinning}; 
        η = 1, τ = 1, 
        estimator = VisitationFrequency(), 
        b = 2)

# Conditional TE with default discretization scheme(s

Calculate transfer entropy from `source` to `response` conditioned on
`cond` using the provided 
`estimator` on a rectangular discretization of a `k + l + m + n`-dimensional 
delay embedding of the input data, using an embedding delay of `τ` across 
all embedding components. `η` is the prediction lag. 

## Arguments
- **`source`**: The source data series.
- **`target`**: The target data series.
- **`cond`**: The conditional data series.
- **`k`**: The dimension of the ``T_{f}`` component of the embedding. 
- **`l`**: The dimension of the ``T_{pp}`` component of the embedding. 
- **`m`**: The dimension of the ``S_{pp}`` component of the embedding. 
- **`n`**: The dimension of the ``C_{pp}`` component of the embedding. 
- **`binning_scheme`**: The binning scheme(s) used to construct the partitions
    over which TE is computed. Must be either one or several instances of 
    `RectangularBinning`s (provided as a vector). TE is computed for each 
    of the resulting partitions.

## Keyword arguments

- **`τ`**: The embedding lag. Default is `τ = 1`.
- **`η`**: The prediction lag. Default is `η = 1`.
- **`estimator`**: The transfer entropy estimator to use. The default 
    is `VisitationFrequency()`.
- **`b`**: Base of the logarithm. The default (`b = 2`) gives the TE in bits.

## More about the embedding

To compute transfer entropy, we need an appropriate delay embedding 
of `source` (``S``), `target` (``T``) and `cond` (``C``). For convenience, define 


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
,``S_{pp}`` denotes the `m`-dimensional set of vectors furnishing the past and present of ``S``, 
and ``C_{pp}`` denotes the `n`-dimensional set of vectors furnishing the past and present of ``C``.
``\\eta`` is the prediction lag. This convenience function uses ``\\tau_1`` = `τ`, 
``\\tau_2`` = `2*τ`, ``\\tau_3`` = `3*τ`, and so on.


Combined, we get the generalised embedding ``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)}, C_{pp}^{(n)})``, 
which is discretized as described below.

## More about discretization


The discretization scheme must be either a single `RectangularBinning` instance, or a vector of 
`RectangularBinning` instances. Run `?RectangularBinning` after loading `CausalityTools` for 
details.


## Transfer entropy computation

TE is then computed as 

```math
\\begin{align}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}, C_{pp}) \\log_{b}{\\left(\\frac{P(T_f | T_{pp}, S_{pp}, C_{pp})}{P(T_f | T_{pp}, C_{pp})}\\right)}
\\end{align}
```
using the provided `estimator` (default = `VisitationFrequency`) for 
each of the discretizations. A vector of the TE estimates for each discretization 
is returned.
"""
function te_cond(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        cond::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int, n::Int,
        binning_scheme::Vector{RectangularBinning}; 
        η = 1, τ = 1, 
        estimator = VisitationFrequency(), 
        b = 2)

    k + l + m + n >= 4 || throw(ArgumentError("`dim = k + l + m + n` must be 4 or higher for conditional TE"))

    pts, vars = te_embed(source, response, cond, k, l, m, n, η = η, τ = τ)

    # Compute TE over the partitions constructed from the provided binning schemes
    # ====================================
    tes = [transferentropy(pts, vars, bs, estimator, b = b) for bs in binning_scheme]

    return tes
end


function te_cond(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        cond::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int, n::Int,
        binning_scheme::RectangularBinning; 
        η = 1, τ = 1, 
        estimator = VisitationFrequency(), 
        b = 2)

    k + l + m + n >= 4 || throw(ArgumentError("`dim = k + l + m + n` must be 4 or higher for conditional TE"))

    pts, vars = te_embed(source, response, cond, k, l, m, n, η = η, τ = τ)

    # Compute TE over the partitions constructed from the provided binning schemes
    # ===========================================================================
    te = transferentropy(pts, vars, binning_scheme, estimator, b = b)

    return te
end

export te_cond