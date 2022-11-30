using Neighborhood

export GaoKannanOhViswanath

"""
    GaoKannanOhViswanath <: MutualInformationEstimator
    GaoKannanOhViswanath(; k = 1, w = 0)

The `GaoKannanOhViswanath` (Shannon) estimator is designed for estimating
mutual information between variables that may be either discrete, continuous or
a mixture of both (Gao et al., 2017).

!!! note "Explicitly convert your discrete data to floats"
    Our implementation of the algorithm demands that all input data are
    numeric and float-valued. If you have discrete data, encode it using
    integers and convert those integers to floats prior to using
    [`mutualinfo`](@ref).

We here extend the estimator to multiple variables (i.e. the multi-information),
which was not treated in Gao et al., (2017).

## Description

The estimator starts by expressing mutual information in terms of the
Radon-Nikodym derivative, and then estimates these derivatives using `k`-nearest neighbor
distances from empirical samples.

The estimator avoids the common issue of having to add noise to data before analysis
due to tied points, which may bias other estimators. Citing their paper, the
estimator *"strongly outperforms natural baselines of discretizing the mixed random
variables (by quantization) or making it continuous by adding a small Gaussian noise."*

!!! note "Resemblance to `KSG1` estimator"
    In Gao et al., (2017), they claim (roughly speaking) that the estimator
    reduces to the [`KSG1`](@ref) estimator for continuous-valued data.
    However, [`KSG1`](@ref) uses the digamma function, while `GaoKannanOhViswanath`
    uses the logarithm instead, so the estimators are not exactly equivalent
    for continuous data.

See also: [`mutualinfo`](@ref).

[^GaoKannanOhViswanath2017]:
    Gao, W., Kannan, S., Oh, S., & Viswanath, P. (2017). Estimating mutual information for
    discrete-continuous mixtures. Advances in neural information processing systems, 30.
"""
Base.@kwdef struct GaoKannanOhViswanath <: MutualInformationEstimator
    k::Int = 1
    w::Int = 0
end

# Note: input datasets must have the same type. Remind the user ot convert in the
# docstring.
function mutualinfo(e::Renyi, est::GaoKannanOhViswanath,
        x, y)
    (; k, w) = est
    joint = Dataset(x, y)
    N = length(joint)
    M = length(x)
    metric = Chebyshev()
    tree_joint = KDTree(joint.data, metric)
    ds_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    tree_x = KDTree(x, metric)
    tree_y = KDTree(y, metric)

    h = 0.0
    for i = 1:N
        # The notation for ρ_{i, xy} in the paper in unclear. They claim in the paper that
        # the estimator reduces to the KSG1 estimator when k̂ == k. Therefore,
        # I assume ρ_{i, xy} is the distance in the *joint* space.
        dmax = ds_joint[i]
        k̂ = dmax == 0 ? inrangecount(tree_joint, joint[i], 0.0) - 1  : k
        h += digamma(k̂) + log(N)

        # They claim in the paper that the estimator reduces to the KSG1 estimator when
        # k̂ == k. However, it only does so when using `digamma`. Their estimator uses
        # `log`, so `GaoKannanOhViswanath` != `KSG1`, but quite close in practice.
        # inrangecount includes the point itself, so we don't need to add 1 inside log
        h -= log(inrangecount(tree_x, x[i], dmax))
        h -= log(inrangecount(tree_y, y[i], dmax))
    end
    h /= N

    return h / log(e.base, ℯ)
end
