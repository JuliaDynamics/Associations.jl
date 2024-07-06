using Neighborhood: Chebyshev, KDTree, Theiler, NeighborNumber
using Neighborhood: bulksearch
using SpecialFunctions: digamma
using DelayEmbeddings: StateSpaceSet, AbstractStateSpaceSet
using Statistics: mean

export KraskovStögbauerGrassberger2, KSG2

"""
    KSG2 <: MutualInformationEstimator
    KraskovStögbauerGrassberger2 <: MutualInformationEstimator
    KraskovStögbauerGrassberger2(; k::Int = 1, w = 0, metric_marginals = Chebyshev())

The `KraskovStögbauerGrassberger2` Shannon mutual information estimator (you can use `KSG2` for
short) is the ``I^{(2)}`` `k`-th nearest neighbor estimator from [Kraskov2004](@cite).

## Usage

- Use with [`association`](@ref) to compute [`MIShannon`](@ref) from input data.

## Keyword arguments

- **`k::Int`**: The number of nearest neighbors to consider. Only information about the
    `k`-th nearest neighbor is actually used.
- **`metric_marginals`**: The distance metric for the marginals for the marginals can be
    any metric from `Distances.jl`. It defaults to `metric_marginals = Chebyshev()`, which
    is the same as in Kraskov et al. (2004).
- **`w::Int`**: The Theiler window, which determines if temporal neighbors are excluded
    during neighbor searches in the joint space. Defaults to `0`, meaning that only the
    point itself is excluded.

## Description

Let the joint StateSpaceSet ``X := \\{\\bf{X}_1, \\bf{X_2}, \\ldots, \\bf{X}_m \\}`` be defined by the
concatenation of the marginal StateSpaceSets ``\\{ \\bf{X}_k \\}_{k=1}^m``, where each ``\\bf{X}_k``
is potentially multivariate. Let ``\\bf{x}_1, \\bf{x}_2, \\ldots, \\bf{x}_N`` be the points
in the joint space ``X``.

The `KraskovStögbauerGrassberger2` estimator first locates, for each ``\\bf{x}_i \\in X``, the
point ``\\bf{n}_i \\in X``, the `k`-th nearest neighbor to ``\\bf{x}_i``, according to the
maximum norm (`Chebyshev` metric). Let ``\\epsilon_i`` be the
distance ``d(\\bf{x}_i, \\bf{n}_i)``.

Consider ``x_i^m \\in \\bf{X}_m``, the ``i``-th point in the marginal space ``\\bf{X}_m``. For each
``\\bf{x}_i^m``, we determine ``\\theta_i^m`` := the number of points ``\\bf{x}_k^m \\in \\bf{X}_m`` that
are a distance less than ``\\epsilon_i`` away from ``\\bf{x}_i^m``. That is, we use the
distance from a query point ``\\bf{x}_i \\in X`` (in the *joint* space) to count neighbors of
``x_i^m \\in \\bf{X}_m`` (in the marginal space).

Shannon mutual information between the variables ``\\bf{X}_1, \\bf{X_2}, \\ldots, \\bf{X}_m`` is
then estimated as

```math
\\hat{I}_{KSG2}(\\bf{X}) =
    \\psi{(k)} -
    \\dfrac{m - 1}{k} +
    (m - 1)\\psi{(N)} -
    \\dfrac{1}{N} \\sum_{i = 1}^N \\sum_{j = 1}^m \\psi{(\\theta_i^j + 1)}
```

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000)
y = rand(rng, 10000)
information(KSG2(; k = 10), x, y) # should be near 0 (and can be negative)
```

## Compatible definitions

- [`MIShannon`](@ref)
"""
struct KraskovStögbauerGrassberger2{M <: MutualInformation, MJ, MM} <: MutualInformationEstimator{M}
    definition::M # the definition of the measure
    k::Int
    w::Int
    metric_joint::MJ # always Chebyshev, otherwise estimator is not valid!
    metric_marginals::MM # can be any metric
end
const KSG2 = KraskovStögbauerGrassberger2

function KraskovStögbauerGrassberger2(definition = MIShannon();
        k::Int = 1,
        w::Int = 0,
        metric_marginals = Chebyshev()
    )
    metric_joint = Chebyshev()
    KraskovStögbauerGrassberger2(definition, k, w, metric_joint, metric_marginals)
end

function information(est::KSG2{<:MIShannon}, x::VectorOrStateSpaceSet...)
    verify_number_of_inputs_vars(est.definition, length(x))

    @assert length(x) >= 2 ||
        error("Need at leats two input StateSpaceSets to compute mutual information between them.")

    (; definition, k, w, metric_joint, metric_marginals) = est
    joint = StateSpaceSet(x...)
    marginals = map(xᵢ -> StateSpaceSet(xᵢ), x)
    M = length(x)
    N = length(joint)

    # `ds[i]` := the distance to k-th nearest neighbor for the point `joint[i]`.
    tree_joint = KDTree(joint, metric_joint)
    knn_idxs, knn_ds = bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))
    idxs = last.(knn_idxs)
    ds = last.(knn_ds)
    ns = [zeros(Int, N) for _ in eachindex(marginals)]
    ϵs = [zeros(N) for _ in eachindex(marginals)]

    for (m, xₘ) in enumerate(marginals)
        eval_dists_to_knns!(ϵs[m], xₘ, idxs, est.metric_marginals)
        marginal_inrangecount!(est, ns[m], xₘ, idxs, ds, m)
    end
    ϵ_maxes = [maximum(x) for x in StateSpaceSet(ϵs...)]
    marginal_nₖs = StateSpaceSet(ns...)

    mi = digamma(k) -
        # The commented-out term appears in Kraskov (2004), but that gives
        # erroneous estimates that do not align with what they show in the
        # paper. Not dividing by k, as done here, fixes the problem,
        # and appears in the KSG2 estimator equation in
        # https://warwick.ac.uk/fac/sci/mathsys/people/students/2019intake/ni/corrected_report_mi_haoranni.pdf,
        # though it is not stated how it is reached.
        (M-1) - #(M - 1) / k  -
        mean(sum(digamma.(nₖ)) for nₖ in marginal_nₖs) +
        (M - 1) * digamma(N)
    return convert_logunit(mi, ℯ, definition.base)
end

function marginal_inrangecount!(est::KraskovStögbauerGrassberger2, ns::Vector{Int},
        xₘ::AbstractStateSpaceSet, knn_idxs, ds, m)
    @assert length(ns) == length(xₘ)
    N = length(xₘ)
    tree = KDTree(xₘ, est.metric_marginals)
    for i = 1:N
        xᵢᵐ = xₘ[i]
        # Add small noise to facilitate ≤ while still using inrangecount
        ϵᵢᵐ = evaluate(est.metric_marginals, xᵢᵐ, xₘ[knn_idxs[i]]) + 1e1*eps()
        # Subtract 1 because `inrangecount` includes the point itself.
        ns[i] = inrangecount(tree, xᵢᵐ, ϵᵢᵐ) - 1
    end
    return ns
end
