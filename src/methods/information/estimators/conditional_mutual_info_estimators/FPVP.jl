
using Neighborhood: bulkisearch, inrangecount
using Neighborhood: Theiler, NeighborNumber, KDTree, Chebyshev
using SpecialFunctions: digamma

export FPVP

"""
    FPVP <: ConditionalMutualInformationEstimator
    FPVP(definition = CMIShannon(); k = 1, w = 0)

The Frenzel-Pompe-Vejmelka-Paluš (or `FPVP` for short)
[`ConditionalMutualInformationEstimator`](@ref) is used to estimate the
conditional mutual information using a `k`-th nearest neighbor approach that is
analogous to that of the [`KraskovStögbauerGrassberger1`](@ref) mutual information
estimator from [Frenzel2007](@citet) and [Vejmelka2008](@citet).

`k` is the number of nearest neighbors. `w` is the Theiler window, which controls the
number of temporal neighbors that are excluded during neighbor searches.

## Usage

- Use with [`association`](@ref) to compute [`ConditionalMutualInformation`](@ref) measure
    from input data.

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000)
y = rand(rng, 10000) .+ x
z = rand(rng, 10000) .+ y
association(FPVP(; k = 10), x, z, y) # should be near 0 (and can be negative)
```

## Compatible definitions

- [`CMIShannon`](@ref)
"""
struct FPVP{M <: ConditionalMutualInformation, MJ, MM} <: ConditionalMutualInformationEstimator{M}
    definition::M
    k::Int
    w::Int
    metric_joint::MJ
    metric_marginals::MM
end

function FPVP(definition = CMIShannon(); k = 1, w = 0)
    # Metrics shouldn't be modified by the user.
    metric_joint = Chebyshev()
    metric_marginals = Chebyshev()
    return FPVP(definition, k, w, metric_joint, metric_marginals)
end

function association(est::FPVP{<:CMIShannon}, x, y, z)
    (; definition, k, w, metric_joint, metric_marginals) = est

    # Ensures that vector-valued inputs are converted to StateSpaceSets, so that
    # building the marginal/joint spaces and neighbor searches are fast.
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    Z = StateSpaceSet(z)
    @assert length(X) == length(Y) == length(Z)
    N = length(X)
    joint = StateSpaceSet(X, Y, Z)
    XZ = StateSpaceSet(X, Z)
    YZ = StateSpaceSet(Y, Z)

    tree_joint = KDTree(joint, metric_joint)
    ds_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    tree_xz = KDTree(XZ, metric_marginals)
    tree_yz = KDTree(YZ, metric_marginals)
    tree_z = KDTree(Z, metric_marginals)

    condmi = 0.0
    for (i, dᵢ) in enumerate(ds_joint)
        # Usually, we subtract 1 because inrangecount includes the point itself,
        # but we'll have to add it again inside the digamma, so just skip it.
        condmi += digamma(k)
        condmi -= digamma(inrangecount(tree_xz, XZ[i], dᵢ))
        condmi -= digamma(inrangecount(tree_yz, YZ[i], dᵢ))
        condmi += digamma(inrangecount(tree_z, Z[i], dᵢ))
    end
    # The "unit" is nats.
    condmi /= N

    return _convert_logunit(condmi, ℯ, definition.base)
end
