export GaussianMI
using StateSpaceSets: Dataset
using StateSpaceSets: dimension, standardize
using LinearAlgebra: eigvals

"""
    GaussianMI <: MutualInformationEstimator
    GaussianMI(; normalize::Bool = false)

`GaussianMI` is a parametric estimator for Shannon mutual information.

## Description

Given ``d_x``-dimensional and ``d_y``-dimensional input data `X` and `Y`,
`GaussianMI` first constructs the ``d_x + d_y``-dimensional joint [`Dataset`](@ref) `XY`.
If `normalize == true`, then we follow the approach in Vejmelka & Palus
(2008)[^Vejmelka2008] and transform each column in `XY` to have zero mean and unit
standard deviation. If `normalize == false`, then the algorithm proceeds without
normalization.

Next, the `C_{XY}`, the correlation matrix for the (normalized) joint data `XY` is
computed. The mutual information estimate

`GaussianMI` assumes the input variables are distributed according to normal
distributions with zero means and unit standard deviations.
Therefore, given ``d_x``-dimensional and ``d_y``-dimensional input data `X` and `Y`,
`GaussianMI` first constructs the joint [`Dataset`](@ref) `XY`, then transforms each
column in `XY` to have zero mean and unit standard deviation, and finally computes
the `\\Sigma`, the correlation matrix for `XY`.

The mutual information estimated (for `normalize == false`) is then estimated as

```math
\\hat{I}^S_{Gaussian}(X; Y) = \\dfrac{1}{2}
\\dfrac{ \\det(\\Sigma_X) \\det(\\Sigma_Y)) }{\\det(\\Sigma))},
```

where we ``\\Sigma_X`` and ``\\Sigma_Y`` appear in ``\\Sigma`` as

```math
\\Sigma = \\begin{matrix}
\\Sigma_{X} & \\Sigma^{'}\\\\
\\Sigma^{'} & \\Sigma_{Y}
\\end{matrix}.
```

If `normalize == true`, then the mutual information is estimated as

```math
\\hat{I}^S_{Gaussian}(X; Y) = -\\dfrac{1}{2} \\sum{i = 1}^{d_x + d_y} \\sigma_i,
```

where ``\\sigma_i`` are the eigenvalues for ``\\Sigma``.

[^Vejmelka2008]:
    Vejmelka, M., & Paluš, M. (2008). Inferring the directionality of coupling with
    conditional mutual information. Physical Review E, 77(2), 026214.
"""
Base.@kwdef struct GaussianMI <: MutualInformationEstimator
    normalize::Bool = false
end

function estimate(measure::MIShannon, est::GaussianMI, x, y)
    X = Dataset(x)
    Y = Dataset(y)
    DX = dimension(X)
    DY = dimension(Y)

    XY = Dataset(X, Y)
    if est.normalize
        Σ = fastcor(standardize(XY))
        σ = eigvals(Σ)
        return -0.5 * sum(log(σᵢ) for σᵢ in σ) / log(measure.e.base, ℯ)
    else
        Σ = fastcor(XY)
        Σx = Σ[1:DX, 1:DX]
        Σy = Σ[DX+1:end, DX+1:end]
        mi = 0.5 * log((det(Σx) * det(Σy)) / det(Σ))
        return mi
    end
end
