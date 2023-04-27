export CorrTest
export CorrTestResult

import HypothesisTests: pvalue

# Note: HypothesisTests already defines CorrelationTest.
# Until https://github.com/JuliaStats/HypothesisTests.jl/issues/290 is addressed and
# HypothesisTests.jl extends pvalue from StatsAPI.jl, we need to use a different name,
# so for now, use `CorrTest`.
"""
    CorrTest <: IndependenceTest
    CorrTest()

An independence test based correlation (for two variables) and partial
correlation (for three variables), as described in Schmidt et al. (2018)[^Schmidt2018].
Uses [`PearsonCorrelation`](@ref) and [`PartialCorrelation`](@ref) internally.

Assumes that the input data are (multivariate) normally distributed. Then
`ρ(X, Y) = 0` implies `X ⫫ Y ` and `ρ(X, Y | 𝐙) = 0` implies `X ⫫ Y | 𝐙`.

## Description

The null hypothesis is `H₀ := ρ(X, Y | 𝐙) = 0`. We use
the approach in Levy & Narula (1978)[^Levy1978] and compute the Z-transformation
of the observed (partial) correlation coefficient ``\\hat{\\rho}_{XY|\\bf{Z}}``:

```math
Z(\\hat{\\rho}_{XY|\\bf{Z}}) =
\\log\\dfrac{1 + \\hat{\\rho}_{XY|\\bf{Z}}}{1 - \\hat{\\rho}_{XY|\\bf{Z}}}.
```

To test the null hypothesis against the alternative hypothesis
`H₁ := ρ(X, Y | 𝐙) > 0`, calculate

```math
\\hat{Z} = \\dfrac{1}{2}\\dfrac{Z(\\hat{\\rho}_{XY|\\bf{Z}}) - Z(0)}{\\sqrt{1/(n - d - 3)}},
```

and compute the two-sided p-value (Schmidt et al., 2018)

```math
p(X, Y | \\bf{Z}) = 2(1 - \\phi(\\sqrt{n - d - 3}Z(\\hat{\\rho}_{XY|\\bf{Z}}))),
```

where ``d`` is the dimension of ``\\bf{Z}`` and ``n`` is the number of samples.
For the pairwise case, the procedure is identical, but set ``\\bf{Z} = \\emptyset``.

## Examples

- [`CorrTest`for independence between normally distributed data](@ref examples_corrtest).

[^Levy1978]:
    Levy, K. J., & Narula, S. C. (1978). Testing hypotheses concerning partial
    correlations: Some methods and discussion. International Statistical Review/Revue Internationale de Statistique, 215-218.

[^Schmidt2018]:
    Schmidt, C., Huegle, J., & Uflacker, M. (2018, July). Order-independent
    constraint-based causal structure learning for gaussian distribution models using
    gpus. In Proceedings of the 30th International Conference on Scientific and
    Statistical Database Management (pp. 1-10).
"""
Base.@kwdef struct CorrTest{M} <: IndependenceTest{M}
    measure::M = nothing
end

"""
    CorrTestResult(pvalue, ρ, z)

A simple struct that holds the results of a [`CorrTest`](@ref) test: the (partial)
correlation coefficient `ρ`, Fisher's `z`, and `pvalue` - the two-sided
p-value for the test.
"""
struct CorrTestResult{R, Z, P}
    ρ::R
    z::Z
    pvalue::P
end
pvalue(r::CorrTestResult) = r.pvalue

function Base.show(io::IO, test::CorrTestResult)
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
        `CorrTest` independence test
        ----------------------------------------------------------------------------------
        H₀: "The first two variables are independent (given the 3rd variable, if relevant)"
        H₁: "The first two variables are dependent (given the 3rd variable, if relevant)"
        ----------------------------------------------------------------------------------
        ρ, (partial) correlation: $(test.ρ)
        p-value (two-sided):      $(test.pvalue)
          $α005
          $α001
          $α0001\
        """
        )
end

const VectorOr1D{D} = Union{AbstractVector, AbstractDataset{D}} where D
function independence(test::CorrTest, x::VectorOr1D, y::VectorOr1D, z::ArrayOrStateSpaceSet...)
    if isempty(z)
        ρ = estimate(PearsonCorrelation(), x, y)
        z = fishers_z(ρ)
        pval = pvalue(test, z, 0, length(x))
        return CorrTestResult(ρ, z, pval)
    else
        # This is already done in `estimate(::PartialCorrelation, ...)`, so may seem
        # redundant, but we need the dimension of Z, so some code duplication occurs here.
        X, Y, Z = construct_partialcor_datasets(x, y, z...)
        D = StateSpaceSet(X, Y, Z)
        cov = fastcov(D)
        precision_matrix = invert_cov(cov)
        ρ = partial_correlation_from_precision(precision_matrix, 1, 2)
        z = fishers_z(ρ)
        pval = pvalue(test, z, dimension(Z), length(x))
        return CorrTestResult(ρ, z, pval)
    end
end

"""
    pvalue(test::CorrTest, z, c::Int, n::Int)

Compute the two-sided p-value for the test of the partial correlation coefficient `p̂`
being zero, where `c` is the cardinality of the conditioning set and `n` is the number of
samples.
"""
function pvalue(test::CorrTest, z, c::Int, n::Int)
    N = Normal(0, 1)
    x = sqrt(n - c - 3) * abs(z)
    return 2*(1 - cdf(N, x))
end
