export GaussianMI
using StateSpaceSets: StateSpaceSet
using StateSpaceSets: dimension, standardize
using LinearAlgebra: eigvals, det

"""
    GaussianMI <: MutualInformationEstimator
    GaussianMI(; normalize::Bool = false)

`GaussianMI` is a parametric estimator for Shannon mutual information.

## Compatible definitions

- [`MIShannon`](@ref)

## Usage

- Use with [`association`](@ref) to compute Shannon mutual information from input data.
- Use with some [`IndependenceTest`](@ref) to test for independence between variables.

## Description

Given ``d_x``-dimensional and ``d_y``-dimensional input data `X` and `Y`,
`GaussianMI` first constructs the ``d_x + d_y``-dimensional joint [`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet) `XY`.
If `normalize == true`, then we follow the approach in Vejmelka & Palus
(2008)[Vejmelka2008](@cite) and transform each column in `XY` to have zero mean and unit
standard deviation. If `normalize == false`, then the algorithm proceeds without
normalization.

Next, the `C_{XY}`, the correlation matrix for the (normalized) joint data `XY` is
computed. The mutual information estimate `GaussianMI` assumes the input variables are distributed according to normal
distributions with zero means and unit standard deviations.
Therefore, given ``d_x``-dimensional and ``d_y``-dimensional input data `X` and `Y`,
`GaussianMI` first constructs the joint [`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet) `XY`, then transforms each
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

## Example 

```julia
using Associations
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000); y = rand(rng, 10000)
association(GaussianMI(), x, y) # should be near 0 (and can be negative)
```
"""
struct GaussianMI{M<:MutualInformation} <: MutualInformationEstimator{M}
    definition::M
    normalize::Bool
end

function GaussianMI(definition=MIShannon(); normalize=true)
    return GaussianMI(definition, normalize)
end

function association(est::GaussianMI{<:MIShannon}, x, y)
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    DX = dimension(X)
    DY = dimension(Y)

    XY = StateSpaceSet(X, Y)
    if est.normalize
        Σ = cor(standardize(XY))
        σ = eigvals(Σ)
        mi = -0.5 * sum(log(σᵢ) for σᵢ in σ)
    else
        Σ = cor(XY)
        Σx = Σ[1:DX, 1:DX]
        Σy = Σ[DX+1:end, DX+1:end]
        mi = 0.5 * log((det(Σx) * det(Σy)) / det(Σ))
    end

    return convert_logunit(mi, ℯ, est.definition.base)
end
