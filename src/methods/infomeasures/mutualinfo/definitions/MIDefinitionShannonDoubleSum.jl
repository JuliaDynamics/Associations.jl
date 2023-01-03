export MIDefinitionShannonDoubleSum

"""
    MIDefinitionShannonDoubleSum <: MutualInformationDefinition

The double-sum definition of discrete Shannon mutual information (see [`MIShannon`](@ref)).

## Description

Assume we observe samples
``\\bar{\\bf{X}}_{1:N_y} = \\{\\bar{\\bf{X}}_1, \\ldots, \\bar{\\bf{X}}_n \\}`` and
``\\bar{\\bf{Y}}_{1:N_x} = \\{\\bar{\\bf{Y}}_1, \\ldots, \\bar{\\bf{Y}}_n \\}`` from
two discrete random variables ``X`` and ``Y`` with finite supports
``\\mathcal{X} = \\{ x_1, x_2, \\ldots, x_{M_x} \\}`` and
``\\mathcal{Y} = y_1, y_2, \\ldots, x_{M_y}``.

In practice, we estimate the double sum `` I(X; Y) = \\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y) log \\left( \\dfrac{p(x, y)}{p(x)p(y)} \\right)``
by replacing the true probabilities ``p(x_i)``, ``p(y_j)``, and ``p(x_i, y_j)`` by their
maximum likelihood empirical estimates ``\\hat{p}(x_i) = \\frac{n(x_i)}{N_x}``,
``\\hat{p}(y_i) = \\frac{n(y_j)}{N_y}``,
and ``\\hat{p}(x_i, x_j) = \\frac{n(x_i)}{N}``, where ``N = N_x N_y``.

Essentially, this is just computing a 2D histogram over the observations (i.e. a
["two-way contingency table"](https://en.wikipedia.org/wiki/Frequency_(statistics)#:~:text=Joint%20frequency-,distributions,-%5Bedit%5D)).
From the (literally) margins of this table/matrix, we can read off the marginal probabilities
and evaluate Shannon mutual information directly by the sum:

```math
\\begin{align*}
\\hat{I}_{DS}(X; Y) &=
 \\sum_{x_i \\in \\mathcal{X}, y_i \\in \\mathcal{Y}} p(x_i, y_j) \\log \\left( \\dfrac{p(x_i, y_i)}{p(x_i)p(y_j)} \\right) \\\\
&= p(x_1, y_1) \\log p(x_1, y_1) + p(x_1, y_2) \\log p(x_1, y_2) + \\cdots + p(x_1, y_{M_y}) \\log p(x_1, y_{M_y}) \\\\
&+ p(x_2, y_1) \\log p(x_2, y_1) + p(x_2, y_2) \\log p(x_2, y_2) + \\cdots + p(x_2, y_{M_y}) \\log p(x_2, y_{M_y}) \\\\
&+ \\cdots \\\\
&+ p(x_{M_x}, y_1) \\log p(x_{M_x}, y_1) + p(x_{M_x}, y_2) \\log p(x_{M_x}, y_2) + \\cdots + p(x_{M_x}, y_1) \\log p(x_{M_x}, y_{M_y}) \\\\
\\end{align*}
```

## Compatibility (discrete estimators)

We currently support the following estimators.

- [`SymbolicPermutation`](@ref). Works with *timeseries* input. Each timeseries
    is separately [`encode`](@ref)d according to its ordinal pattern, and probabilities
    are estimated of relative frequencies of the ordinal patterns.
- [`Dispersion`](@ref). Works with *timeseries* input. Each timeseries
    is separately [`encode`](@ref)d according to its ordinal pattern, and probabilities
    are estimated of relative frequencies of the ordinal patterns.

Although accurate for discrete data, this method can be a bit slow to compute in practice,
because to estimate the joint probabilities ``p(x_i, y_j)``, we have to explicitly
compute the joint distribution histogram, which can be quite expensive. An alternative,
although less accurate, is the 3-entropies formulation [`MIDefinitionShannonH3`](@ref).

!!! info "Potential for more implementations!"
    The double-sum definition is in principle compatible with any
    [`ProbabilitiesEstimator`](@ref) for which the [`outcome_space`](@ref) is known,
    given the input data. However, because mutual information is a function of multiple
    inputs, this requires specialized implementations that encode input data into their
    outcome representations, so that we can explicitly compute the contingency table.

    A plethora of implementations are possible. Pull requests are welcome!
"""
struct MIDefinitionShannonDoubleSum <: MutualInformationDefinition end
