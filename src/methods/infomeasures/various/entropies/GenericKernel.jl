export GenericKernel

"""
    GenericKernel <: DifferentialInfoEstimator
    GenericKernel(bandwidth = Silverman(), kernel::MultivariateKernel = NormalIsotropic())

A generic, multivariate plug-in estimator for entropies based on kernel density estimation
(KDE) that can in principle be used to compute any differential entropy.

Data should be standardized to zero mean and unit variance before applying `GenericKernel`.

The `bandwidth` may be set manually, or to one of the rule-of-thumbs listed below.

## Description

Assume we have samples ``\\{x_1, x_2, \\ldots, x_N \\}`` from a continuous random variable
``X \\in \\mathbb{R}^d`` with support ``\\mathcal{X}`` and density function
``f : \\mathbb{R}^d \\to \\mathbb{R}``.

`GenericKernel` estimates, for each ``x_i`` in the sample, the point-wise densities
``\\hat{f}(x_i)`` using the given `kernel` and `bandwidth`, then computes a resubstitution
estimate for the entropy. We support the following resubstitution estimates.

### [Shannon](@ref) differential entropy

```math
H(X) = \\int_{\\mathcal{X}} f(x) \\log f(x) dx = \\mathbb{E}[-\\log(f(X))]
```

is estimated by replacing the expectation with the sample average ([^Diks2017])

```math
\\hat{H}(X) = -\\dfrac{1}{N}\\sum_{i = 1}^N \\log \\hat{f}(x).
```

## Compatible kernels

- [`NormalIsotropic`](@ref).
- [`Epanechnikov`](@ref)

## Bandwidth rule-of-thumbs

- [`Silverman`](@ref)
- [`DiksFang`](@ref)

[^Diks2017].
    Diks, C., & Fang, H. (2017). Transfer entropy for nonparametric granger causality
    detection: an evaluation of different resampling methods. Entropy, 19(7), 372.
"""
struct GenericKernel{K, B} <: DifferentialInfoEstimator
    bandwidth::B
    kernel::K
    function GenericKernel(
            bandwidth::B = DiksFang(4.8),
            kernel::K = NormalIsotropic()) where {B, K}
        new{K, B}(bandwidth, kernel)
    end
end
bandwidth(r::Real, x::AbstractStateSpaceSet) = r # convenience for manual settings

function entropy(e::Renyi, est::GenericKernel, x::AbstractStateSpaceSet)
    bw = bandwidth(est.bandwidth, x)
    ρs = densities_at_points(est.kernel, x, bw)
    e.q ≈ 1 || error("Renyi entropy with q = $(e.q) not implemented for `GenericKernel`")
    return sum(log0.(e.base, ρs)) / length(x)
end
