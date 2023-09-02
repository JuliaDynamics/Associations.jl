export GaussianMI
using StateSpaceSets: StateSpaceSet
using StateSpaceSets: dimension, standardize
using LinearAlgebra: eigvals

"""
    GaussianMI <: MutualInformationEstimator
    GaussianMI(; normalize::Bool = false)

`GaussianMI` is a parametric estimator for Shannon mutual information.

## Description

Given ``d_x``-dimensional and ``d_y``-dimensional input data `X` and `Y`,
`GaussianMI` first constructs the ``d_x + d_y``-dimensional joint [`StateSpaceSet`](@ref) `XY`.
If `normalize == true`, then we follow the approach in Vejmelka & Palus
(2008)[Vejmelka2008](@cite) and transform each column in `XY` to have zero mean and unit
standard deviation. If `normalize == false`, then the algorithm proceeds without
normalization.

Next, the `C_{XY}`, the correlation matrix for the (normalized) joint data `XY` is
computed. The mutual information estimate `GaussianMI` assumes the input variables are distributed according to normal
distributions with zero means and unit standard deviations.
Therefore, given ``d_x``-dimensional and ``d_y``-dimensional input data `X` and `Y`,
`GaussianMI` first constructs the joint [`StateSpaceSet`](@ref) `XY`, then transforms each
column in `XY` to have zero mean and unit standard deviation, and finally computes
the `\\Sigma`, the correlation matrix for `XY`.

The mutual information estimated (for `normalize == false`) is then estimated as

```math
\\hat{I}^S_{Gaussian}(X; Y) = \\dfrac{1}{2}
\\dfrac{ \\det(\\Sigma_X) \\det(\\Sigma_Y)) }{\\det(\\Sigma))},
```

where we ``\\Sigma_X`` and ``\\Sigma_Y`` appear in ``\\Sigma`` as

```math
\\Sigma = \\begin{bmatrix}
\\Sigma_{X} & \\Sigma^{'}\\\\
\\Sigma^{'} & \\Sigma_{Y}
\\end{bmatrix}.
```

If `normalize == true`, then the mutual information is estimated as

```math
\\hat{I}^S_{Gaussian}(X; Y) = -\\dfrac{1}{2} \\sum_{i = 1}^{d_x + d_y} \\sigma_i,
```

where ``\\sigma_i`` are the eigenvalues for ``\\Sigma``.
"""
Base.@kwdef struct GaussianMI <: MutualInformationEstimator
    normalize::Bool = false
end

function estimate(measure::MIShannon, est::GaussianMI, x, y)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    DX = dimension(X)
    DY = dimension(Y)

    XY = StateSpaceSet(X, Y)
    if est.normalize
        Σ = fastcor(standardize(XY))
        σ = eigvals(Σ)
        mi = -0.5 * sum(log(σᵢ) for σᵢ in σ)
    else
        Σ = fastcor(XY)
        Σx = Σ[1:DX, 1:DX]
        Σy = Σ[DX+1:end, DX+1:end]
        mi = 0.5 * log((det(Σx) * det(Σy)) / det(Σ))
    end

    return convert_logunit(mi, ℯ, measure.e.base)
end
