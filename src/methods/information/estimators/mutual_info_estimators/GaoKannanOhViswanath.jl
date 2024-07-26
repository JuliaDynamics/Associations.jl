using Neighborhood

export GaoKannanOhViswanath

"""
    GaoKannanOhViswanath <: MutualInformationEstimator
    GaoKannanOhViswanath(; k = 1, w = 0)

The `GaoKannanOhViswanath` (Shannon) estimator is designed for estimating
Shannon mutual information between variables that may be either discrete, continuous or
a mixture of both [GaoKannanOhViswanath2017](@cite).

## Compatible definitions

- [`MIShannon`](@ref)

## Usage

- Use with [`association`](@ref) to compute Shannon mutual information from input data.
- Use with some [`IndependenceTest`](@ref) to test for independence between variables.

## Description

The estimator starts by expressing mutual information in terms of the
Radon-Nikodym derivative, and then estimates these derivatives using `k`-nearest neighbor
distances from empirical samples.

The estimator avoids the common issue of having to add noise to data before analysis
due to tied points, which may bias other estimators. Citing their paper, the
estimator *"strongly outperforms natural baselines of discretizing the mixed random
variables (by quantization) or making it continuous by adding a small Gaussian noise."*

!!! warn "Implementation note"
    In [GaoKannanOhViswanath2017](@citet), they claim (roughly speaking) that the estimator
    reduces to the [`KraskovStögbauerGrassberger1`](@ref) estimator for continuous-valued data.
    However, [`KraskovStögbauerGrassberger1`](@ref) uses the digamma function, while `GaoKannanOhViswanath`
    uses the logarithm instead, so the estimators are not exactly equivalent
    for continuous data.

    Moreover, in their algorithm 1, it is clearly not the case that the method falls
    back on the `KraskovStögbauerGrassberger1` approach. The `KraskovStögbauerGrassberger1` estimator uses `k`-th neighbor distances in
    the *joint* space, while the `GaoKannanOhViswanath` algorithm selects the maximum
    `k`-th nearest distances among the two marginal spaces, which are in general not the
    same as the `k`-th neighbor distance in the joint space (unless both marginals are
    univariate). Therefore, our implementation here differs slightly from algorithm 1 in
    `GaoKannanOhViswanath`. We have modified it in a way that mimics [`KraskovStögbauerGrassberger1`](@ref) for
    continous data. Note that because of using the `log` function instead of `digamma`,
    there will be slight differences between the methods. See the source code for more
    details.

!!! note "Explicitly convert your discrete data to floats"
    Even though the `GaoKannanOhViswanath` estimator is designed to handle discrete data,
    our implementation demands that all input data are `StateSpaceSet`s whose data points
    are floats. If you have discrete data, such as strings or symbols, encode them using
    integers and convert those integers to floats before passing them to
    [`association`](@ref).

## Examples

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000); y = rand(rng, 10000)
association(GaoKannanOhViswanath(; k = 10), x, y) # should be near 0 (and can be negative)
```
"""
struct GaoKannanOhViswanath{M <: MutualInformation} <: MutualInformationEstimator{M}
    definition::M
    k::Int
    w::Int 
end

function GaoKannanOhViswanath(definition = MIShannon(); k = 1, w = 0)
    return GaoKannanOhViswanath(definition, k, w)
end
# TODO: We here extend the estimator to multiple variables (i.e. the multi-information),
# which was not treated in Gao et al., (2017).

# Note: input StateSpaceSets must have the same type. Remind the user ot convert in the
# docstring.
function association(est::GaoKannanOhViswanath{<:MIShannon}, x, y)
    (; definition, k, w) = est
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    joint = StateSpaceSet(X, Y)

    N = length(joint)
    M = length(x)
    metric = Chebyshev()
    tree_joint = KDTree(joint, metric)
    tree_x = KDTree(X, metric)
    tree_y = KDTree(Y, metric)
    ds_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    ds_x = last.(bulksearch(tree_x, X, NeighborNumber(k), Theiler(w))[2])
    ds_y = last.(bulksearch(tree_y, Y, NeighborNumber(k), Theiler(w))[2])

    mi = 0.0
    for i = 1:N
        # The notation for ρ_{i, xy} in the paper in unclear. They claim in the paper that
        # the estimator reduces to the KSG1 estimator when k̂ == k, i.e. when Therefore,
        # I assume they mean ρ_{i, xy} is the distance in the *joint* space. However, then
        # the estimator fails for continuous data.
        dmax = max(ds_x[i], ds_y[i])
        if dmax == 0
            # Don't subtract 1 here for inrangecount. We want precisely the number
            # if points such that the distance is zero.
            k̂ = inrangecount(tree_joint, joint[i], 0.0)
        else
            dmax = ds_joint[i]
            k̂ = k
        end
        mi += digamma(k̂) + log(N)

        # They claim in the paper that the estimator reduces to the KSG1 estimator when
        # k̂ == k. However, it only does so when using `digamma`. Their estimator uses
        # `log`, so `GaoKannanOhViswanath` != `KraskovStögbauerGrassberger1`, but quite close in practice.
        # inrangecount includes the point itself, so we don't need to add 1 inside log
        nx = inrangecount(tree_x, X[i], dmax) - 1
        ny = inrangecount(tree_y, Y[i], dmax) - 1
        mi -= log(nx)
        mi -= log(ny)
    end
    # The "unit" is nats.
    mi /= N
    return _convert_logunit(mi, ℯ, definition.base)
end
