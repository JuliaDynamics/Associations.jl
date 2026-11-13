export GaussianMI
using StateSpaceSets: StateSpaceSet
using StateSpaceSets: dimension, standardize
using LinearAlgebra: det

"""
    GaussianMI <: MutualInformationEstimator
    GaussianMI(definition = MIShannon(); normalize::Bool = true)

`GaussianMI` is a parametric estimator for Shannon mutual information.

## Compatible definitions

- [`MIShannon`](@ref)

## Usage

- Use with [`association`](@ref) to compute Shannon mutual information from input data.
- Use with some [`IndependenceTest`](@ref) to test for independence between variables.

## Description

Given ``d_x``-dimensional and ``d_y``-dimensional input data `X` and `Y`, `GaussianMI`
models the joint data `XY` as Gaussian. It computes ``\\Sigma``, the correlation matrix of
`XY`, with blocks

```math
\\Sigma = \\begin{bmatrix}
\\Sigma_{X} & \\Sigma^{'}\\\\
\\Sigma^{'T} & \\Sigma_{Y}
\\end{bmatrix},
```

and estimates the mutual information as

```math
\\hat{I}^S_{Gaussian}(X; Y) = \\dfrac{1}{2}
\\log \\dfrac{ \\det(\\Sigma_X) \\det(\\Sigma_Y) }{\\det(\\Sigma)}.
```

If `normalize == true`, then each column of `XY` is first transformed to have zero mean
and unit standard deviation, following Vejmelka & Palus (2008)[Vejmelka2008](@cite).
Since ``\\Sigma`` is a correlation matrix, this does not change the estimate.

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
    # `Σ` is a correlation matrix, so standardizing the columns first does not change it.
    Σ = cor(est.normalize ? standardize(XY) : XY)
    Σx = Σ[1:DX, 1:DX]
    Σy = Σ[DX+1:end, DX+1:end]
    mi = 0.5 * log((det(Σx) * det(Σy)) / det(Σ))

    return convert_logunit(mi, ℯ, est.definition.base)
end
