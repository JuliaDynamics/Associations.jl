using Random

export PC

"""
    PC <: GraphAlgorithm
    PC(unconditional_test, conditional_test; α = 0.05, max_depth = nothing,
        maxiters_orient = ∞)

The PC algorithm (Spirtes et al., 2000)[^Spirtes2000], implemented as described in
Kalisch & Bühlmann (2008)[^Kalisch2008].

## Arguments

- **`unconditional_test`**: An [`IndependenceTest`](@ref) that uses a pairwise,
    nondirectional [`AssociationMeasure`](@ref) measure, e.g. [`CorrTest`](@ref)
    or [`SurrogateTest`](@ref) with the [`MIShannon`](@ref) measure.
- **`conditional_test`**: An [`IndependenceTest`](@ref) that uses a conditional,
    nondirectional [`AssociationMeasure`](@ref), , e.g. [`CorrTest`](@ref)
    or [`SurrogateTest`](@ref) with the [`CMIShannon`](@ref) measure.

## Keyword arguments

- **`α::Real ∈ (0, 1)`**. The significance level of the test.
- **`max_depth::Union{Int, Nothing} = nothing`**. The significance level of the test.
- **`maxiters_orient::Real`**. The maximum number of times to apply the orientation
    rules. By default, there is not limit (i.e. `maxiters_orient = ∞`)

!!! info "Directional measures will not give meaningful answers"
    During the skeleton search phase, if a significance association between two nodes are
    is found, then a bidirectional edge is drawn between them. The generic implementation of
    `PC` therefore doesn't currently handle directional measures such as
    [`TEShannon`](@ref). The reason is that if a  directional relationship `X → Y` exists
    between two nodes `X` and `Y`, then the algorithm would first draw a bidirectional
    arrow between `X` and `Y` when analysing the direction
    `X → Y`, and then removing it again when analysing in the direction `Y → X` (a similar
    situation would also occur for the conditional stage). This will be fixed in a
    future release. For now, use nondirectional measures, e.g. [`MIShannon`](@ref) and
    [`CMIShannon`](@ref)!

## Description

When used with [`infer_graph`](@ref) on some input data `x`, the `PC` algorithm
performs the following steps:

1. Initialize an empty fully connected graph `g` with `N` nodes, where `N` is the number
    of variables and `x[i]` is the data for the `i`-th node.
2. Reduce the fully connected `g` to a skeleton graph by performing pairwise
    [`independence`](@ref) tests between all vertices using `unconditional_test`. Remove
    any edges where adjacent vertices are found to be independent according to the test
    (i.e. the null hypothesis of independence cannot be rejected at significance level
    `1 - α`).
3. Thin the skeleton `g` by conditional [`independence`](@ref) testing. If
    `x[i] ⫫ x[j] | x[Z]` for some set of variables `Z` (not including `i` and `j`)
    according to `conditional_test` (i.e. the null hypothesis of conditional independence
    cannot be rejected at significance level `1 - α`), then the edge between `i` and `j` is
    removed, and we record the separating set S(i, j) = Z. Independence tests are first
    performed for conditioning sets of size 1, and repeated for conditioning sets of
    increasing size, which in most cases limits the number of tests needed.  The separating
    sets `S(i, j)`, which records which variables were in the conditioning set that
    rendered variables `i` and `j` independent, are recorded.
    If `max_depth` is an integer, then this procedure is performed on conditioning
    sets of sizes `1:max_depth`, and if `max_depth == nothing`, then all possible
    conditioning set sizes are potentially used.
4. Create a directed graph `dg` from `g` by replacing every
    undirected edge `X - Y` in `g` by the bidirectional edge `X ↔ Y` (i.e.
    construct two directional edges `X → Y` and `Y → X`). Orientiation rules 0-3
    are then repeatedly applied to `dg` until no more edges can be oriented:
    - Rule 0 (orients v-structures): `X ↔ Y ↔ Z` becomes `X → Y ← Z` if `Y` is not in the
        separating set `S(X, Z)`.
    - Rule 1 (prevents new v-structures): `X → Y ↔ Z` becomes `X → Y → Z` if `X` and `Z`
        are not adjacent.
    - Rule 2 (avoids cycles): `X → Y → Z ↔ X` becomes `X → Y → Z ← X`.
    - Rule 3: To avoid creating cycles or new v-structures, whenever `X - Y → Z`,
        `X - W → Z`, and `X - Z` but there is no edge between `Y` and `W`, turn the
        undirected `X - Z` edge into the directed edge `X → Z`.

The resulting directed graph (a `SimpleDiGraph` from
[Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl)) is then returned.


## Examples

- [Inferring a causal graph from a chain of logistic maps](@ref pc_robust_example)

!!! info
    The "PC" algorithm is named after the *first names* of the authors, **P**eter Spirtes
    and **C**lark Glymour.

[^Kalisch2008]:
    Kalisch, M., & Bühlmann, P. (2008). Robustification of the PC-algorithm for directed
    acyclic graphs. Journal of Computational and Graphical Statistics, 17(4), 773-789.
[^Spirtes2000]:
    Spirtes, P., Glymour, C. N., Scheines, R., & Heckerman, D. (2000). Causation, prediction,
    and search. MIT press.
"""
struct PC{U, C, A, N} <: GraphAlgorithm
    unconditional_test::U
    conditional_test::C
    α::A
    maxdepth::N
    maxiters_orient::Int

    function PC(
            unconditional_test::IndependenceTest,
            conditional_test::IndependenceTest;
            α::A = 0.05, maxdepth::N = nothing,
            maxiters_orient::Int = 500) where {A, N <: Union{Int, Nothing}}
        0 < α < 1 || throw(ArgumentError("α must be on `(0, 1)`. α = 0.05 is commonly used"))
        U = typeof(unconditional_test)
        C = typeof(conditional_test)
        new{U, C, A, N}(unconditional_test, conditional_test, α, maxdepth, maxiters_orient)
    end
end

include("skeleton.jl")
#include("skeleton_directed.jl")
include("cpdag.jl")

# TODO: this is only designed for non-directed measures. Use type system?
function infer_graph(algorithm::PC, x; verbose = false)
    sg, separating_sets = skeleton(algorithm, x; verbose)
    dg = cpdag(algorithm, sg, separating_sets; verbose)

    return sg, dg
end
