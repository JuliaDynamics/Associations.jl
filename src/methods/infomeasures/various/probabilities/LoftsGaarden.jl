export LoftsGaarden
import ComplexityMeasures: outcome_space

"""
    Loftsgaarden <: ProbabilitiesEstimator

The `Loftsgaarden` probabilities estimator is based on the `k`-th nearest neighbor
density estimatio from Loftsgaarden & Quesenberry (1965).

It estimates probabilities by first estimating the density locally at each sample
point `xᵢ` using the distance from `xᵢ` to its `k`-th nearest neighbor. The density
distribution over the sample points is then normalized to form probabilities.

## Outcome space

The outcome space `Ω` for `LoftsGaarden` is the indices of the input data, `1:length(x)`.
The reason to not return the data points themselves is because duplicate data points may
not have same probabilities (due to having different neighbors).

[^Loftsgaarden1965]:
    Loftsgaarden, D. O., & Quesenberry, C. P. (1965). A nonparametric estimate of a
    multivariate density function. The Annals of Mathematical Statistics, 36(3), 1049-1051.
"""
Base.@kwdef struct LoftsGaarden{M} <: ProbabilitiesEstimator
    k::Int = 5
    w::Int = 0
    metric::M = Euclidean()
end

function probabilities_and_outcomes(est::LoftsGaarden, x::AbstractDataset{D}) where D
    Probabilities(point_densities(est, x)), 1:length(x)
end

outcome_space(x::AbstractDataset, ::LoftsGaarden) = 1:length(x)

function point_densities(est::LoftsGaarden, x::AbstractDataset{D}) where D
    (; k, w, metric) = est
    N = length(x)
    bᵥ = ComplexityMeasures.ball_volume(D)

    # The bandwidth `bws[i]` for the point `x[i]` is the distance to the `k`-th nearest
    # neighbor of `x[i]`. The local density around, in contrast, in formed by the `kmax`
    # nearest neighbors.
    tree = KDTree(x, metric)
    ds = last.(bulksearch(tree, x, NeighborNumber(k), Theiler(w))[2])

    densities = zeros(N)
    for (i, dᵢ) in enumerate(ds)
        densities[i] = point_density(est, dᵢ, N, bᵥ)
    end

    return densities
end

point_density(est::LoftsGaarden, dᵢ, N, bᵥ) = est.k / (N*bᵥ*dᵢ)
