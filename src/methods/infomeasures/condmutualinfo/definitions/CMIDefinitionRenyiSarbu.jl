export CMIRenyiSarbu

"""
    CMIRenyiSarbu <: ConditionalMutualInformationDefinition
    CMIRenyiSarbu()

A definition of Rényi conditional mutual information (Sarbu, 2014[^Sarbu2014]).

## Description

Assume we observe three discrete random variables ``X``, ``Y`` and ``Z``.
Sarbu (2014) defines discrete conditional Rényi mutual information as the conditional
Rényi ``\\alpha``-divergence between the conditional joint probability mass function
``p(x, y | z)`` and the product of the conditional marginals, ``p(x |z) \\cdot p(y|z)``:

```math
I(X, Y; Z)^R_q =
\\dfrac{1}{q-1} \\sum_{z \\in Z} p(Z = z)
\\log \\left(
    \\sum{x \\in X}\\sum{y \\in Y}
    \\dfrac{p(x, y|z)^q}{\\left( p(x|z)\\cdot p(y|z) \\right)^{q-1}}
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
struct CMIRenyiSarbu <: ConditionalMutualInformationDefinition end
