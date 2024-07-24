export PartialCorrelation

"""
    PartialCorrelation <: AssociationMeasure

The correlation of two variables, with the effect of a set of conditioning
variables removed.

## Usage

- Use with [`association`](@ref) to compute the raw partial correlation coefficient.
- Use with [`independence`](@ref) to perform a formal hypothesis test for
    correlated-based conditional independence.

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
struct PartialCorrelation <: CorrelationMeasure end

min_inputs_vars(::PartialCorrelation) = 3
max_inputs_vars(::PartialCorrelation) = Inf

# Compatibility with `independence`
function association(::PartialCorrelation, x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet,
        conds::ArrayOrStateSpaceSet...)
    X, Y, Z = construct_partialcor_datasets(x, y, conds...)
    D = StateSpaceSet(X, Y, Z)
    cov_matrix = cov(D)
    precision_matrix = invert_cov(cov_matrix)
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

function invert_cov(cov_matrix::AbstractMatrix)
    if det(cov_matrix) ≈ 0.0
        # If the determinant of the covariance matrix is zero, then the
        # Moore-Penrose pseudo-inverse is used.
        rtol = sqrt(eps(real(float(one(eltype(cov_matrix))))))
        return pinv(cov_matrix; rtol)
    else
        return inv(cov_matrix)
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
