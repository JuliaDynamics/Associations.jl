using Distances
"""
    te(driver, response;
        method = :transferoperator_grid,
        E = nothing, v = nothing,
        ϵ = nothing, n_ϵ = 5, lag = 1, dim = 3, η = 1,
        k1 = 2, k2 = 3, metric = Chebyshev(),
        which_is_surr = -1, surr_type = aaft)

Compute transfer entropy from a `driver` time series to a `response` time
series.

## Embedding instructions
Specify embedding and instructions for how to compute marginal
entropies are with `E::AbstractEmbedding` and `v::TEVars`. If either
is missing, the data is embedded using the provided forward `lag`, embedding
dimension `dim` and backward lag `η`.

`method` sets the transfer entropy algorithm to use. There are two types of
estimators: grid-based approaches, and nearest-neighbor based approaches.

## Grid-based estimators

- `:transferoperator_grid`. Grid-based transfer entropy estimator based on the
    transfer operator.
- `:visitfreq`. Simple visitation frequency based estimator.

For all grid based algorithms, `ϵ` sets the binsize. If `ϵ = nothing`, then the
algorithm uses a bin size corresponding to a number of subdivisions `N`
along each axis so that `N ≦ npoints^(1/(dim+1))`. Transfer entropy is then
computed as an average over `n_ϵ` different bin sizes corresponding to
`N-1` to `N` subdivisions along each axis.


## Nearest neighbor based estimator.
- `:kraskov`: A nearest-neighbor based estimator, which computes transfer
    entropy as the sum of two mutual information terms. `k1` and `k2` sets the
    number of nearest neighbors for those terms. `metric` sets the distance
    metric (must be valid metric from `Distances.jl`).

## Surrogate analysis
A surrogate analysis will be run if `which_is_surr` is set to 0, 1 or 2.
- `which_is_surr = 0` will construct a surrogate for both time series.
- `which_is_surr = 1` will construct a surrogate for the driver time series.
- `which_is_surr = 2` will construct a surrogate for the response time series.

"""

function te(driver, response;
        method = :transferoperator_grid,
        E = nothing, v = nothing,
        ϵ = nothing, n_ϵ = 5, lag = 1, dim = 3, η = 1,
        k1 = 2, k2 = 3, metric = Chebyshev(),
        which_is_surr = -1, surr_type = aaft)

    valid_methods = [:transferoperator_grid,
                    :kraskov,
                    :visitfreq]

    if !(method ∈ valid_methods)
        error("Transfer entropy method $method not valid.")
    end

    if which_is_surr == 1
        driver = surr_type(driver)
    elseif which_is_surr == 2
        response = surr_type(response)
    elseif which_is_surr == 0
        driver = surr_type(driver)
        response = surr_type(response)
    elseif which_is_surr >= 3
        error("which_is_surr $which_is_surr is not valid.")
    end

    # Compute transfer entropy analogously to how Krakovskà et al. (2018)
    # computes conditional mutual information.
    if !(typeof(E) == AbstractEmbedding) || !(typeof(v) == TEVars)
        if dim == 3
            E = embed([driver, response], [2, 2, 1], [lag, 0, 0])
            v = TEVars([1], [2], [3])
        elseif dim == 4
            E = embed([driver, response], [2, 2, 2, 1], [lag, 0, -η, 0])
            v = TEVars([1], [2, 3], [4])
        elseif dim == 5
            E = embed([driver, response], [2, 2, 2, 2, 1], [lag, 0, -η, -2*η, 0])
            v = TEVars([1], [2, 3, 4], [5])
        elseif dim == 6
            E = embed([driver, response], [2, 2, 2, 2, 2, 1], [lag, 0, -η, -2*η, -3*η, 0])
            v = TEVars([1], [2, 3, 4, 5], [6])
        elseif dim == 7
            E = embed([driver, response], [2, 2, 2, 2, 2, 2, 1], [lag, 0, -η, -2*η, -3*η, -4*η, 0])
            v = TEVars([1], [2, 3, 4, 5, 6], [7])
        elseif dim == 8
            E = embed([driver, response], [2, 2, 2, 2, 2, 2, 2, 1], [lag, 0, -η, -2*η, -3*η, -4*η, -5*η, 0])
            v = TEVars([1], [2, 3, 4, 5, 6, 7], [8])
        end
    end

    if !(ϵ == nothing)
        ϵs = [ϵ for i = 1:dim]
    else
        ϵs = zeros(Float64, n_ϵ, dim)
        # We a range of bin sizes along each axis.
        for d in 1:dim
            dmin = minimum(E.points[d, :])
            dmax = maximum(E.points[d, :])
            maxrange = dmax - dmin
            max_number_of_bins = ceil(Int, npoints(E)^(1/(dim+1)))
            max_bin_size = maxrange/max_number_of_bins
            min_bin_size = maxrange/(max_number_of_bins - 1)
            ϵs[:, d] = logspace(log(10, min_bin_size), log(10, max_bin_size), n_ϵ)
        end
    end

    tes = zeros(Float64, n_ϵ)
    for i = 1:n_ϵ
        if method == :transferoperator_grid
            tes[i] = transferentropy_transferoperator_visitfreq(E, ϵs[i, :], v)

        elseif method == :kraskov
            tes[i] = transferentropy_kraskov(E, k1, k2, v, metric = metric)

        elseif method == :visitfreq
            tes[i] = transferentropy_visitfreq(E, ϵs[i, :], v)

        elseif method == :transferoperator_triang
            warn("Invariantizing embedding before calculating transfer matrix...")
            t = triangulate(invariantize(E))
            TO = transferoperator_approx(t, n_pts = n_randpts, sample_randomly = sample_randomly)
            invdist = left_eigenvector(TO)
            tes[i] = transferentropy_transferoperator_triang(t, invdist, n_ϵ, n_reps)
        end
    end
    # If we only use one binsize, we can't take the integral.
    if !(ϵ == nothing) || method == :kraskov
        return tes[1]
    else
        ∫te(1:n_ϵ, tes)
    end
end
