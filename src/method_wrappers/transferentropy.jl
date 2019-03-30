
"""
    transferentropy(driver, response;
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
function transferentropy(driver::AbstractVector{T}, response::AbstractVector{T};
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
            tes[i] = transferentropy(E, v, RectangularBinning(ϵs[i]), TransferOperatorGrid())

        elseif estimator ∈ (:transferentropy_kraskov,
                            :tekraskov,
                            :tekNN)
            tes[i] = transferentropy_kraskov(E, k1, k2, v, metric = distance_metric)

        elseif estimator ∈ (:transferentropy_visitfreq,
                            :tefreq)
            #tes[i] = transferentropy_visitfreq(E, ϵs[i], v)
            tes[i] = transferentropy(E, v, RectangularBinning(ϵs[i]), VisitationFreq())

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
te = transferentropy
export transferentropy, te
