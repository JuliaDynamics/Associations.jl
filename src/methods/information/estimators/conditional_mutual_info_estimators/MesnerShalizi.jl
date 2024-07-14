export MesnerShalizi
export MesnerShalisi
"""
    MesnerShalizi <: ConditionalMutualInformationEstimator
    MesnerShalizi(definition = CMIShannon(); k = 1, w = 0)

The `MesnerShalizi` [`ConditionalMutualInformationEstimator`](@ref) is designed for
data that can be mixtures of discrete and continuous data [Mesner2020](@cite).

`k` is the number of nearest neighbors. `w` is the Theiler window, which controls the
number of temporal neighbors that are excluded during neighbor searches.

## Usage

- Use with [`association`](@ref) to compute [`CMIShannon`](@ref) from input data.

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000)
y = rand(rng, 10000) .+ x
z = rand(rng, 10000) .+ y
association(MesnerShalizi(; k = 10), x, z, y) # should be near 0 (and can be negative)
```

## Compatible definitions

- [`CMIShannon`](@ref)
"""
struct MesnerShalizi{M <: ConditionalMutualInformation, ME} <: ConditionalMutualInformationEstimator{M}
    definition::M
    k::Int
    w::Int
    metric::ME
end
function MesnerShalizi(definition = CMIShannon(); k = 1, w = 0)
    # Metric shouldn't be modified by the user.
    metric = Chebyshev()
    return MesnerShalisi(definition, k, w, metric)
end

function MesnerShalisi(args...; kwargs...)
    # Silently deprecate.
    return MesnerShalizi(args...; kwargs...)
end

function association(est::MesnerShalizi{<:CMIShannon}, x, y, z)
 
    (; definition, k, w, metric) = est
    X = StateSpaceSet(x)
    Y = StateSpaceSet(y)
    Z = StateSpaceSet(z)
    joint = StateSpaceSet(X, Y, Z)
    XZ = StateSpaceSet(X, Z)
    YZ = StateSpaceSet(Y, Z)
    N = length(joint)
    M = 3
    tree_joint = KDTree(joint, metric)
    ds_joint = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])
    tree_xz = KDTree(XZ, metric)
    tree_yz = KDTree(YZ, metric)
    tree_z = KDTree(Z, metric)

    condmi = 0.0
    for i = 1:N
        # The notation for ρ_{i, xy} in the paper in unclear. They claim in the paper that
        # the estimator reduces to the KSG1 estimator when k̂ == k. Therefore,
        # I assume ρ_{i, xy} is the distance in the *joint* space.
        # TODO: this might not be correct..."
        dmax = ds_joint[i]
        k̂ = dmax == 0 ? inrangecount(tree_joint, joint[i], 0.0) - 1  : k
        condmi += digamma(k̂)
        # Simulate ≤ by adding smallest possible nudge.
        condmi -= log(inrangecount(tree_xz, XZ[i], dmax + eps()))
        condmi -= log(inrangecount(tree_yz, YZ[i], dmax + eps()))
        condmi += log(inrangecount(tree_z, Z[i], dmax + eps()))
    end
    # The "unit" is nats.
    condmi /= N
    return _convert_logunit(condmi, ℯ, definition.base)
end
