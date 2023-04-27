export partial_correlation
export PartialCorrelation

"""
    PartialCorrelation <: AssociationMeasure

The correlation of two variables, with the effect of a set of conditioning
variables removed.

## Usage

- Use with [`independence`](@ref) to perform a formal hypothesis test for
    conditional dependence.
- Use with [`partial_correlation`](@ref) to compute the raw correlation coefficient.

## Description

There are several ways of estimating the partial correlation. We follow the
[matrix inversion method](https://en.wikipedia.org/wiki/Partial_correlation), because
for [`StateSpaceSet`](@ref)s, we can very efficiently compute the required
joint covariance matrix ``\\Sigma`` for the random variables.

Formally, let ``X_1, X_2, \\ldots, X_n`` be a set of ``n`` real-valued random variables.
Consider the joint precision matrix,``P = (p_{ij}) = \\Sigma^-1``. The partial
correlation of any pair of variables ``(X_i, X_j)``, given the remaining variables
``\\bf{Z} = \\{X_k\\}_{i=1, i \\neq i, j}^n``, is defined as

```math
\\rho_{X_i X_j | \\bf{Z}} = -\\dfrac{p_ij}{\\sqrt{ p_{ii} p_{jj} }}
```

In practice, we compute the estimate

```math
\\hat{\\rho}_{X_i X_j | \\bf{Z}} =
-\\dfrac{\\hat{p}_ij}{\\sqrt{ \\hat{p}_{ii} \\hat{p}_{jj} }},
```

where ``\\hat{P} = \\hat{\\Sigma}^{-1}`` is the sample precision matrix.
"""
struct PartialCorrelation <: AssociationMeasure end

"""
    partial_correlation(x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet,
        z::VectorOrStateSpaceSet...)

Compute the [`PartialCorrelation`](@ref) between `x` and `y`, given `z`.
"""
function partial_correlation(x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet, z::ArrayOrStateSpaceSet...)
    return estimate(PartialCorrelation(), x, y, z...)
end

# Compatibility with `independence`
function estimate(::PartialCorrelation, x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet,
        conds::ArrayOrStateSpaceSet...)
    X, Y, Z = construct_partialcor_datasets(x, y, conds...)
    D = StateSpaceSet(X, Y, Z)
    cov = fastcov(D)
    precision_matrix = invert_cov(cov)
    return partial_correlation_from_precision(precision_matrix, 1, 2)
end

function construct_partialcor_datasets(x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet,
        conds::ArrayOrStateSpaceSet...)
    dimension(x) == 1 || throw(ArgumentError("Input `x` must be 1-dimensional"))
    dimension(y) == 1 || throw(ArgumentError("Input `y` must be 1-dimensional"))
    X, Y = StateSpaceSet(x), StateSpaceSet(y)
    Z = StateSpaceSet(conds...)
    return X, Y, Z
end

function estimate(measure::PartialCorrelation, est::Nothing, x, y, z)
    return estimate(measure, x, y, z)
end


function invert_cov(cov::AbstractMatrix)
    if det(cov) ≈ 0.0
        # If the determinant of the covariance matrix is zero, then the
        # Moore-Penrose pseudo-inverse is used.
        rtol = sqrt(eps(real(float(one(eltype(cov))))))
        return pinv(cov; rtol)
    else
        return inv(cov)
    end
end

"""
    partial_correlation_from_precision(P, i::Int, j::Int)

Given a precision matrix `P`, compute the partial correlation of variables `i` and `j`
conditional on all other variables.
"""
function partial_correlation_from_precision(P::AbstractMatrix, i::Int, j::Int)
    return -P[i, j] / sqrt(P[i, i]*P[j, j])
end
