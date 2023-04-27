"""
    LoftsgaardenH <: ProbabilitiesEstimator

The `LoftsGaardenH` Shannon entropy estimator is based on the `k`-th nearest neighbor
density estimation from Loftsgaarden & Quesenberry (1965).

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

function entropy(e::Renyi, est::LoftsGaarden, x)
    ρs = point_densities(est, StateSpaceSet(x))
end
