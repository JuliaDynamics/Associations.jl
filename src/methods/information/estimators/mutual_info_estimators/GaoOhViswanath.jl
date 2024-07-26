using ComplexityMeasures: ball_volume
export GaoOhViswanath

"""
    GaoOhViswanath <: MutualInformationEstimator

The `GaoOhViswanath` is a mutual information estimator based on nearest neighbors,
and is also called the bias-improved-KSG estimator, or BI-KSG, by [Gao2018](@cite).

## Compatible definitions

- [`MIShannon`](@ref)

## Usage

- Use with [`association`](@ref) to compute Shannon mutual information from input data.
- Use with some [`IndependenceTest`](@ref) to test for independence between variables.

## Description

The estimator is given by

```math
\\begin{align*}
\\hat{H}_{GAO}(X, Y)
&= \\hat{H}_{KSG}(X) + \\hat{H}_{KSG}(Y) - \\hat{H}_{KZL}(X, Y) \\\\
&= \\psi{(k)} +
    \\log{(N)} +
    \\log{
        \\left(
            \\dfrac{c_{d_{x}, 2} c_{d_{y}, 2}}{c_{d_{x} + d_{y}, 2}}
        \\right)
    } - \\\\
    & \\dfrac{1}{N} \\sum_{i=1}^N \\left( \\log{(n_{x, i, 2})} + \\log{(n_{y, i, 2})} \\right)
\\end{align*},
```

where ``c_{d, 2} = \\dfrac{\\pi^{\\frac{d}{2}}}{\\Gamma{(\\dfrac{d}{2} + 1)}}`` is the
volume of a ``d``-dimensional unit ``\\mathcal{l}_2``-ball.

## Example 

```julia
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 10000); y = rand(rng, 10000)
association(GaoOhViswanath(; k = 10), x, y) # should be near 0 (and can be negative)
```
"""
struct GaoOhViswanath{M <: MutualInformation, MJ, MM} <: MutualInformationEstimator{M}
    definition::M
    k::Int
    w::Int
    metric_joint::MJ
    metric_marginals::MM
end

function GaoOhViswanath(definition = MIShannon(); 
        k = 1, w = 0, 
        metric_joint = Euclidean(),
        metric_marginals = Euclidean()
    )
    return GaoOhViswanath(definition, k, w, metric_joint, metric_marginals)
end

function association(est::GaoOhViswanath{<:MIShannon}, x::VectorOrStateSpaceSet...)
    verify_number_of_inputs_vars(est.definition, length(x))

    @assert length(x) >= 2 ||
        error("Need at leats two input StateSpaceSets to compute mutual information between them.")
    (; k, w, metric_joint, metric_marginals) = est
    joint = StateSpaceSet(x...)
    marginals = map(xᵢ -> StateSpaceSet(xᵢ), x)
    M = length(x)
    N = length(joint)

    # `ds[i]` := the distance to k-th nearest neighbor for the point `joint[i]`.
    # In the paper, for method 1, ϵᵢ = 2*ds[i].
    tree_joint = KDTree(joint, metric_joint)
    ds = last.(bulksearch(tree_joint, joint, NeighborNumber(k), Theiler(w))[2])

    # `marginals_nₖs[i][k]` contains the number of points within radius `ds[i]` of
    # the point `marginals[k][i]`.
    ns = [zeros(Int, N) for m in eachindex(marginals)]
    s = 0.0
    for (m, xₘ) in enumerate(marginals)
        marginal_inrangecount!(est, ns[m], xₘ,  ds)
    end
    marginal_nₖs = StateSpaceSet(ns...)

    bvₘs = prod(ball_volume(dimension(xₘ)) for xₘ in marginals)

    mi = digamma(k) +
        (M - 1) * log(N) +
        log(bvₘs / ball_volume(dimension(joint))) -
        (1 / N) * sum(sum(log.(nₖ)) for nₖ in marginal_nₖs)
    return _convert_logunit(mi, ℯ, est.definition.base)
end

function marginal_inrangecount!(est::GaoOhViswanath, ns, xₘ, ds)
    @assert length(ns) == length(xₘ)
    tree = KDTree(xₘ, est.metric_marginals)
    @inbounds for i in eachindex(xₘ)
        # Subtract 1 because `inrangecount` includes the point itself.
        ns[i] = inrangecount(tree, xₘ[i], ds[i]) - 1
    end
    return ns
end

# ## KSG

# ```math
# \\begin{align}
# \\hat{H}_{KSG}(X, Y)
# &= \\hat{H}_{KSG}(X) + \\hat{H}_{KSG}(Y) - \\hat{H}_{KZL}(X, Y) \\\\
# \\end{align}
# ```

# where ``\\hat{H}_{KSG}`` is the Kraskov-Stögbauer-Grassberger (KSG) marginal entropy
# estimator estimator

# ```math
# \\hat{H}_{KSG} = \\dfrac{1}{N} \\sum_{i}^N \\log{
#     \\left(-\\psi{(n_{x, i, \\infty} + 1)} +
#     \\psi{(N)} +
#     \\log{(c_{d_{x}, \\infty})} -
#     d_x \\log{(\\rho_{k, i, \\infty})}
#     \\right)}
# ```

# and ``\\hat{H}_{KZL}`` is the Kozachenko-Leonenko (KZL) marginal entropy estimator

# ```math
# \\hat{H}_{KZL} = \\dfrac{1}{N} \\sum_{i}^N \\log{
#     \\left(
#         \\dfrac{N \\cdot c_{d, p} \\cdot (\\rho_{k, i, p})^d}{k}
#     \\right)} + \\log{(k)} - \\psi{(k)},
# ```

# and ``\\psi{(x)}`` is the digamma function.
