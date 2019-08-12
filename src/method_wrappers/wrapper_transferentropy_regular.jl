import StatsBase.mean

"""
    standard_te(driver, response;
        estimator = :tetogrid,
        E = nothing, v = nothing,
        ϵ = nothing, n_ϵ = 5, η = 1, τ = 1, dim = 3,
        k1 = 2, k2 = 3, distance_metric = Chebyshev(),
        which_is_surr = :none,
        surr_func = aaft,
        min_numbins = nothing,
        max_numbins = nothing)

Compute transfer entropy from a `driver` time series to a `response` time
series.


## Method
`method` sets the transfer entropy estimator to use. There are two types of
estimators: grid-based approaches, and nearest-neighbor based approaches.


## Embedding instructions
Specify embedding and instructions for how to compute marginal
entropies are with `E::AbstractEmbedding` and `v::TEVars`. If either
is missing, the data is embedded using the provided forward prediction lag `η`,
with embedding dimension `dim` and embedding lag `η`.

## Grid-based estimators

- `:transferentropy_transferoperator_grid`, or `:tetogrid`. Grid-based transfer entropy
    estimator based on the transfer operator.
- `::transferentropy_visitfreq`, or `:tefreq`. Simple visitation frequency based estimator.

For the grid based TE estimators, `ϵ` sets the bin size. If `ϵ == nothing`, then
the algorithm uses a bin size corresponding to a number of subdivisions `N`
along each axis so that `N ≦ npoints^(1/(dim+1))`. Transfer entropy is then
computed as an average over `n_ϵ` different bin sizes corresponding to
`N-1` to `N` subdivisions along each axis. If `ϵ` is specified and consists
of multiple bin sizes, an average
TE estimate over the bin sizes is returened. If `ϵ` is a single value,
then TE is estimated for that bin size only.


## Nearest neighbor based estimator.
- `:transferentropy_kraskov`, `:tekraskov`, or `:tekNN`. A nearest-neighbor
    based estimator that computes transfer entropy as the sum of two mutual
    information terms. `k1` and `k2` sets the
    number of nearest neighbors used to estimate the mutual information terms.
    `distance_metric` sets the distance metric (must be an instance
    of a valid metric from `Distances.jl`).

## Surrogate analysis
A surrogate analysis will be run if `which_is_surr` is set to
    `:both`, `:driver`, `:response`. Default is `:none`.
- `which_is_surr = :both` will replace both time series with a surrogate.
- `which_is_surr = :driver` will replace the driver time series with a surrogate.
- `which_is_surr = :response` will replace the response time series with a surrogate.

The type of surrogate must be either `randomshuffle`, `randomphases`,
`randomamplitudes`, `aaft` or `iaaft`.
"""
function standard_te(driver::AbstractVector{T}, response::AbstractVector{T};
        dim = 3,
        τ = 1,
        η = 1,
        estimator = :tetogrid,
        E = nothing, v = nothing,
        ϵ = nothing, n_ϵ = 5,
        k1 = 2, k2 = 3, distance_metric = Chebyshev(),
        which_is_surr = :none,
        surr_func = aaft,
        nbin_range = 2,
        min_numbins = nothing,
        max_numbins = nothing) where T
 
    estims_te_grid = [:transferentropy_transferoperator_grid, :tetogrid]
    estims_te_kNN = [:transferentropy_kraskov, :tekraskov, :tekNN]
    estims_te_freq = [:transferentropy_visitfreq, :tefreq]
    valid_estimators = estims_te_grid ∪ estims_te_kNN ∪ estims_te_freq

    if dim < 3
        error("Dimension $dim is too low. Must be at least 3.")
    end

    if !(estimator ∈ valid_estimators)
        error("Transfer entropy method $estimator not valid.")
    end

    valid_surr_funcs = [TimeseriesSurrogates.randomshuffle,
                        TimeseriesSurrogates.randomphases,
                        TimeseriesSurrogates.randomamplitudes,
                        TimeseriesSurrogates.aaft,
                        TimeseriesSurrogates.iaaft]

    if !(surr_func ∈ valid_surr_funcs)
        error("Surrogate type $surr_func is not valid.")
    end
    if which_is_surr == :driver
        driver[:] = surr_func(driver)
    elseif which_is_surr == :response
        response[:] = surr_func(response)
    elseif which_is_surr == :both
        driver[:] = surr_func(driver)
        response[:] = surr_func(response)
    elseif !(which_is_surr == :none)
        error("$which_is_surr is not a valid surrogate choice." )
    end
    D = Dataset(driver, response)
    # Compute transfer entropy analogously to how Krakovskà et al. (2018)
    # computes conditional mutual information.
    if !(typeof(E) == Embeddings.AbstractEmbedding) || !(typeof(v) == TEVars)
        if dim == 3
            E = customembed(D, Positions(2, 2, 1), Lags(η, 0, 0))
            v = TEVars(Tf = [1], Tpp = [2], Spp = [3])
        elseif dim == 4
            E = customembed(D, Positions(2, 2, 2, 1), Lags(η, 0, -τ, 0))
            v = TEVars(Tf = [1], Tpp = [2, 3], Spp = [4])
        elseif dim == 5
            E = customembed(D, Positions(2, 2, 2, 2, 1), Lags(η, 0, -τ, -2*τ, 0))
            v = TEVars(Tf = [1], Tpp = [2, 3, 4], Spp = [5])
        elseif dim == 6
            E = customembed(D, Positions(2, 2, 2, 2, 2, 1), Lags(η, 0, -τ, -2*τ, -3*τ, 0))
            v = TEVars(Tf = [1], Tpp = [2, 3, 4, 5], Spp = [6])
        elseif dim == 7
            E = customembed(D, Positions(2, 2, 2, 2, 2, 2, 1), Lags(η, 0, -τ, -2*τ, -3*τ, -4*τ, 0))
            v = TEVars(Tf = [1], Tpp = [2, 3, 4, 5, 6], Spp = [7])
        elseif dim == 8
            E = customembed(D, Positions(2, 2, 2, 2, 2, 2, 2, 1), Lags(η, 0, -τ, -2*τ, -3*τ, -4*τ, -5*τ, 0))
            v = TEVars(Tf = [1], Tpp = [2, 3, 4, 5, 6, 7], Spp = [8])
        end
    end

    if ϵ == nothing
        ϵ = zeros(Float64, n_ϵ, dim)
        # We a range of bin sizes along each axis.
        for d in 1:dim
            dmin = Base.minimum(E.reconstructed_pts[:, d])
            dmax = Base.maximum(E.reconstructed_pts[:, d])
            maxrange = dmax - dmin

            # If minimum and maximum number of bins are not provided,
            # use the Palus recommendation of N^(1/dim) as a max.
            if !(typeof(min_numbins) == Int && typeof(max_numbins) == Int)
                max_numbins = ceil(Int, length(E.reconstructed_pts)^(1/(dim+1)))
                min_numbins = max_numbins - nbin_range
            end

            max_bin_size = maxrange/max_numbins
            min_bin_size = maxrange/min_numbins

            ϵ[:, d] = LinRange(log(10, min_bin_size), log(10, max_bin_size), n_ϵ)
        end
        binsizes = [ϵ[i, :] for i = 1:n_ϵ]
        ϵs = binsizes
    else
        if length(ϵ) == 1
            ϵs = [ϵ]
        else
            ϵs = ϵ
        end
    end

    tes = zeros(Float64, length(ϵs))

    # Loop over the bin sizes
    for i = 1:length(ϵs)
        if estimator ∈ (:transferentropy_transferoperator_grid,
                        :tetogrid)
            #tes[i] = transferentropy_transferoperator_grid(E, ϵs[i], v)
            tes[i] = transferentropy(E.reconstructed_pts, v, RectangularBinning(ϵs[i]), TransferOperatorGrid())

        elseif estimator ∈ (:transferentropy_kraskov,
                            :tekraskov,
                            :tekNN)
            tes[i] = transferentropy_kraskov(E, k1, k2, v, metric = distance_metric)

        elseif estimator ∈ (:transferentropy_visitfreq,
                            :tefreq)
            #tes[i] = transferentropy_visitfreq(E, ϵs[i], v)
            tes[i] = transferentropy(E.reconstructed_pts, v, RectangularBinning(ϵs[i]), VisitationFreq())

        elseif estimator == :transferoperator_triang
                warn("high level wrapper for $estimator is not implemented yet.")
        end

    end

    # If we only use one binsize, we can't take the integral. This is always
    # the case for the kNN based estimators, because they don't use a bin size.
    if estimator ∈ estims_te_kNN
        return tes[1]
    else
        return mean(tes)
    end
end

"""
    te_reg(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int; 
        η = 1, τ = 1, 
        estimator = VisitationFrequency(), 
        n_subdivs = 1,
        b = 2)

## TE estimation with default discretization scheme(s)

Calculate transfer entropy from `source` to `response` using the provided 
`estimator` on a rectangular discretization of a `k + l + m`-dimensional 
delay embedding of the input data, using an embedding delay of `τ` across 
all embedding components. `η` is the prediction lag. 

## Arguments
- **`source`**: The source data series.
- **`target`**: The target data series.
- **`k`**: The dimension of the ``T_{f}`` component of the embedding. 
- **`l`**: The dimension of the ``T_{pp}`` component of the embedding. 
- **`m`**: The dimension of the ``S_{pp}`` component of the embedding. 

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
of `source` (``S``) and `target` (``T``). For convenience, define 


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


Combined, we get the generalised embedding ``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)})``, 
which is discretized as described below.

## More about discretization

To compute TE, we coarse-grain the `k+l+m`-dimensional generalised embedding 
``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)})`` into hyperrectangular boxes.
The magnitude of the TE may be biased by the particular choice of binning scheme,
so we compute TE across a number of different box sizes, determined as follows.

Let ``L`` be the number of observations in ``S`` (and ``T``). The coarsest box size is 
determined by subdiving the i-th coordinate axis into 
``N = ceiling(L^\\frac{1}{k + l + m + 1})`` intervals of equal lengths, 
resulting in a box size of ``|max(dim_{i}) - min(dim_{i})|/N``. The next box size 
is given by ``|max(dim_{i}) - min(dim_{i})|/(N+1)``, then 
``|max(dim_{i}) - min(dim_{i})|/(N+2)``, and so on, until the finest box 
size which is given by ``|max(dim_{i}) - min(dim_{i})|/(N+N_{subdivs})``. 

## Transfer entropy computation

TE is then computed as 

```math
\\begin{align}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}) \\log_{b}{\\left(\\frac{P(T_f | T_{pp}, S_{pp})}{P(T_f | T_{pp})}\\right)}
\\end{align}
```
using the provided `estimator` (default = `VisitationFrequency`) for 
each of the discretizations. A vector of the TE estimates for each discretization 
is returned.
"""
function te_reg(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int; 
        η = 1, τ = 1, 
        estimator = VisitationFrequency(), 
        n_subdivs = 2,
        b = 2)

    k + l + m >= 3 || throw(ArgumentError("`dim = k + l + m` must be 3 or higher for regular TE"))

    pts, vars = te_embed(source, response, k, l, m, η = η, τ = τ)

    # Determine appropriate binnings from time series length (roughly according to 
    # Krakovska et al. (2018)'s recommendations)
    L = length(source) 
    
    # one subdivision coarser than Krakovska as the coarest partition
    nbins_coarsest = ceil(Int, L^(1/(k + l + m + 1))) - 1 
    
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
    te_reg(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int,
        binning_scheme::Vector{RectangularBinning}; 
        η = 1, τ = 1, 
        estimator = VisitationFrequency(), 
        b = 2)

# TE with user-provided discretization scheme(s)

Calculate transfer entropy from `source` to `response` using the provided 
`estimator` on discretizations constructed by the provided `binning_scheme`(s) 
over a `k + l + m`-dimensional delay embedding of the input data, 
using an embedding delay of `τ` across all embedding components. 
`η` is the prediction lag. 

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

- **`τ`**: The embedding lag. Default is `τ = 1`.
- **`η`**: The prediction lag. Default is `η = 1`.
- **`estimator`**: The transfer entropy estimator to use. The default 
    is `VisitationFrequency()`.
- **`b`**: Base of the logarithm. The default (`b = 2`) gives the TE in bits.

## More about the embedding

To compute transfer entropy, we need an appropriate delay embedding 
of `source` (``S``) and `target` (``T``). For convenience, define 


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


Combined, we get the generalised embedding ``\\mathbb{E} = (T_f^{(k)}, T_{pp}^{(l)}, S_{pp}^{(m)})``, 
which is discretized as described below.

## More about discretization

The discretization scheme must be either a single `RectangularBinning` instance, or a vector of 
`RectangularBinning` instances. Run `?RectangularBinning` after loading `CausalityTools` for 
details.

## Transfer entropy computation

TE is then computed as 

```math
\\begin{align}
TE_{S \\rightarrow T} = \\int_{\\mathbb{E}} P(T_f, T_{pp}, S_{pp}) \\log_{b}{\\left(\\frac{P(T_f | T_{pp}, S_{pp})}{P(T_f | T_{pp})}\\right)}
\\end{align}
```
using the provided `estimator` (default = `VisitationFrequency`) for 
each of the discretizations. A vector of the TE estimates for each discretization 
is returned.
"""
function te_reg(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int,
        binning_scheme::Vector{RectangularBinning},
        η = 1, τ = 1, 
        estimator::TransferEntropyEstimator = VisitationFrequency(), 
        b = 2)

    k + l + m >= 3 || throw(ArgumentError("`dim = k + l + m` must be 3 or higher for regular TE"))

    pts, vars = te_embed(source, response, k, l, m, η = η, τ = τ)

    # Compute TE over the partitions constructed from the provided binning schemes
    # ====================================
    tes = [transferentropy(pts, vars, bs, estimator, b = b) for bs in binning_scheme]

    return tes
end

function te_reg(source::AbstractArray{<:Real, 1}, 
        response::AbstractArray{<:Real, 1}, 
        k::Int, l::Int, m::Int,
        binning_scheme::RectangularBinning,
        η = 1, τ = 1, 
        estimator::TransferEntropyEstimator = VisitationFrequency(), 
        b = 2)

    k + l + m >= 3 || throw(ArgumentError("`dim = k + l + m` must be 3 or higher for regular TE"))

    pts, vars = te_embed(source, response, k, l, m, η = η, τ = τ)

    # Compute TE over the partitions constructed from the provided binning schemes
    # ====================================
    te = transferentropy(pts, vars, binning_scheme, estimator, b = b)

    return te
end

export te_reg