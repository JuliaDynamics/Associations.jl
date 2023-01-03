export MIDefinitionRenyiSarbu


"""
    MIDefinitionRenyiSarbu <: ConditionalMutualInformationDefinition
    MMIDefinitionRenyiSarbu()

A definition of discrete Rényi mutual information (Sarbu, 2014[^Sarbu2014]).

!!! warn "This is slooow"
   Because this definition doesn't lend itself to the very convenient decomposition
   into marginal entropies, using it is quite slow.

## Description

Assume we observe three discrete random variables ``X``, ``Y`` and ``Z``.
Sarbu (2014) defines discrete Rényi mutual information as the
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y)`` and the product of the conditional marginals, ``p(x) \\cdot p(y)``:

```math
I(X, Y; Z)^R_q =
\\dfrac{1}{q-1}
\\log \\left(
    \\sum{x \\in X}\\sum{y \\in Y}
    \\dfrac{p(x, y)^q}{\\left( p(x)\\cdot p(y) \\right)^{q-1}}
\\right)
```

!!! note "Only valid for fixed outcome spaces"
    This estimator is only well defined for [`ProbabilitiesEstimator`](@ref)s with fixed
    outcomes spaces. This is synonymous with having to estimate probabilities in the
    *joint* space, the marginalizing to obtain joint and marginal probability distributions
    with the same number of elements. This works out-of-the-box for some estimators
    like [`ValueHistogram`](@ref) with [`FixedRectangularBinning`](@ref), while for
    estimators such as [`SymbolicPermutation`](@ref), we provide specalized implementations.

[^Sarbu2014]: Sarbu, S. (2014, May). Rényi information transfer: Partial Rényi transfer
    entropy and partial Rényi mutual information. In 2014 IEEE International Conference
    on Acoustics, Speech and Signal Processing (ICASSP) (pp. 5666-5670). IEEE.
"""
struct MIDefinitionRenyiSarbu <: MutualInformationDefinition end
