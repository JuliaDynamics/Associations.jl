"""
    CorrelationTest <: IndependenceTest

An parametric independence test based on correlation (for two variables) and partial
correlation (for three variables), as described in Schmidt et al. (2018)[^Schmidt2018].

## Assumptions

Assumes that the input data are (multivariate) normally distributed.

[^Schmidt2018]:
    Schmidt, C., Huegle, J., & Uflacker, M. (2018, July). Order-independent
    constraint-based causal structure learning for gaussian distribution models using
    gpus. In Proceedings of the 30th International Conference on Scientific and
    Statistical Database Management (pp. 1-10).
"""
struct CorrelationTest <: IndependenceTest end

"""
    CorrelationTestResult(pvalue, ρ, z)

A simple struct that holds the results of a [`Correlation`](@ref) test: the `pvalue`,
the (partial) correlation coefficient `ρ`, and Fisher's `z`.
"""
struct CorrelationTestResult{R, Z, P}
    pvalue::P
    ρ::R
    z::Z
end

function Base.show(io::IO, test::CorrelationTestResult)
    α005 = pvalue(test) < 0.05 ?
        "α = 0.05:  ✓ Evidence favors dependence" :
        "α = 0.05:  ✖ Independence cannot be rejected"
    α001 = pvalue(test) < 0.01 ?
        "α = 0.01:  ✓ Evidence favors dependence" :
        "α = 0.01:  ✖ Independence cannot be rejected"
    α0001 = pvalue(test) < 0.001 ?
        "α = 0.001: ✓ Evidence favors dependence" :
        "α = 0.001: ✖ Independence cannot be rejected"

    print(io,
        """\
        `Correlation` independence test
        ----------------------------------------------------------------------------------
        H₀: "The first two variables are independent (given the 3rd variable, if relevant)"
        H₁: "The first two variables are dependent (given the 3rd variable, if relevant)"
        ----------------------------------------------------------------------------------
        ρ, (partial) correlation: $(test.ρ)
        p-value:                  $(test.pvalue)
          $α005
          $α001
          $α0001\
        """
        )
end
pvalue(r::CorrelationTestResult) = r.pvalue

"""
    pvalue(z, c::Int, n::Int)

Compute the p-value for the test of the partial correlation coefficient `p̂ᵢⱼ` being zero,
where `c` is the cardinality of the conditioning set and `n` is the number of
samples.
"""
function pvalue(z, c::Int, n::Int)
    N = Normal(0, 1)
    x = sqrt(n - c - 3) * abs(z)
    return 2*(1 - cdf(N, x))
end

"""
    fishers_z(p̂ᵢⱼ)

Compute Fisher's z-transform on the sample partial correlation coefficient `p̂ᵢⱼ` (computed
as the correlation between variables `i` and `j`, given the remaining variables):

```math
Z(V_i, V_j | \\bf{S}) = \\dfrac{1}{2}
\\log{\\left(
    \\dfrac{1 + \\hat{p}(V_i, V_j | \\bf{S})}{1 - \\hat{p}(V_i, V_j | \\bf{S})}
\\right) }
```
"""
fishers_z(p̂ᵢⱼ) = 0.5 * log((1 + p̂ᵢⱼ) / (1 - p̂ᵢⱼ))

function independence(test::Correlation, s, t, conds...)
    S, T, C = Dataset(s), Dataset(t), Dataset(conds...)
    D = Dataset(S, T, C)
    cov = fastcov(D)
    return _independence(test, length(D), cov, dimension(C))
end

function _independence(test::Correlation, N::Int, cov, dim_cond::Int)
    # For computing partial correlations, we follow the matrix inversion
    # approach outlined in https://en.wikipedia.org/wiki/Partial_correlation.
    # If the determinant of the covariance matrix is zero, then the
    # Moore-Penrose pseudo-inverse is used.
    rtol = sqrt(eps(real(float(one(eltype(cov))))))
    precision_matrix = det(cov) ≈ 0.0 ? pinv(cov; rtol) : inv(cov)
    ρ = partialcor(precision_matrix, 1, 2)
    z = fishers_z(ρ)
    pval = pvalue(z, dim_cond, N)
    return CorrelationTestResult(pval, ρ, z)
end
