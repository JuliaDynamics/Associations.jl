
using Graphs: SimpleGraphs, SimpleGraph, SimpleDiGraph, Edge
using Graphs: nv, edges, has_edge, inneighbors

"""
    cpdag(alg::PCRobust, skeleton::SimpleGraph, separating_sets::Dict{Edge, Vector{Int}}) → dg::SimpleDiGraph

Orient edges in the `skeleton` graph using the given `separating_sets` using
algorithm 2 in Kalisch & Bühlmann (2008[^Kalisch2008]) and return the directed graph `cpdag`.

## Description

First, a directed graph `dg` is constructed from `skeleton`, such that every
undirected edge `X - Y` in `skeleton` is replaced by the bidirectional edge `X ↔ Y` in `dg`.
In practices, for each `X ↔ Y`, we construct two directional edges `X → Y` and `Y → X`.

Orientiation rules 0-3 are then applied to `dg`:
- Rule 0 (orients v-structures): `X ↔ Y ↔ Z` becomes `X → Y ← Z` if `Y` is not in the
    separating set `S(X, Z)`.
- Rule 1 (prevents new v-structures): `X → Y ↔ Z` becomes `X → Y → Z` if `X` and `Z`
    are not adjacent.
- Rule 2 (avoids cycles): `X → Y → Z ↔ X` becomes X → Y → Z ← X`
- Rule 3: To avoid creating cycles or new v-structures,
       X                 X
     | | |             | | |
    Y  |  W  becomes  Y  |  W
     ↘ | ↙             ↘ ↓ ↙
       Z                 Z

[^Kalisch2008]:
    Kalisch, M., & Bühlmann, P. (2008). Robustification of the PC-algorithm for directed
    acyclic graphs. Journal of Computational and Graphical Statistics, 17(4), 773-789.
"""
function cpdag(alg::PCRobust, skeleton::SimpleDiGraph,
        separating_sets::Dict{SimpleEdge, Vector{Int}}; verbose = false)
    # Convert the skeleton to a directed graph where (Xᵢ - Xⱼ) becomes (Xᵢ → Xⱼ, Xⱼ → Xᵢ).
    dg = skeleton

    # Orientiation rules are described in a plethora of books and papers in the literature.
    # I found most of them, including the original paper on the PC algorithm, hard to
    # understand, due to either terse and/or ambiguous language. The best
    # non-ambiguous description of the rules I found were in Kalisch & Bühlmann (2008).
    # These are applied here. No effort has been done to make this as efficient as possible.
    verbose && println("Transforming skeleton graph into cpdag...")
    verbose && println("⅂ Level 0:")
    rule0!(alg, dg, separating_sets; verbose)

    imax = 20
    rules = ones(3)
    i = 1
    while any(rules .> 0) && i < imax
        verbose && println("⅂ Level 1:")
        # Each application of a rule returns the number of edges that were removed.
        # If no edge was removed in any step, we're done.
        rules[1] = rule1!(alg, dg; verbose)
        rules[2] = rule2!(alg, dg; verbose)
        rules[3] = rule3!(alg, dg; verbose)
        i += 1
    end
    return dg
end

# Rule 0. Orient v-structures. `X ↔ Y ↔ Z` becomes `X → Y ← Z` if `Y` is not in the
# separating set `S(X, Z)`.
function rule0!(alg::PCRobust, dg::SimpleDiGraph, separating_sets::Dict;
        verbose = false)
    n_removed = 0
    for edge in edges(dg)
        !has_edge(dg, reverse(edge)) && continue # Only consider bidirectional edges
        X, Y = edge.src, edge.dst
        for Z in setdiff(inneighbors(dg, Y), X)
            e = SimpleEdge(X, Z)
            if has_edge(dg, Z, Y) &&
                    nonadjacent(dg, X, Z) &&
                    haskey(separating_sets, e) &&
                    Y ∉ separating_sets[e]
                rem_edge!(dg, Z, X)
                rem_edge!(dg, Z, Y)
                verbose && println("⅂  (Rule 0) Removed $Z → $X and $Z → $Y")
                n_removed += 1
            end
        end
    end
end

function rule1!(alg::PCRobust, dg::SimpleDiGraph; verbose = false)
    n_removed = 0
    for edge in edges(dg)
        has_edge(dg, reverse(edge)) && continue # Only consider unidirectional edges
        X, Y = edge.src, edge.dst
        for Z in setdiff(inneighbors(dg, Y), X)
            if has_edge(dg, Z, Y) && nonadjacent(dg, X, Z)
                rem_edge!(dg, Z, Y)
                verbose && println("⅂  (Rule 1) Removed $Z → $Y ")
                n_removed += 1
            end
        end
    end
    return n_removed
end

function rule2!(alg::PCRobust, dg::SimpleDiGraph; verbose = false)
    n_removed = 0

    for edge in edges(dg)
        Y, Z = edge.src, edge.dst
        # X → Y → Z must be a chain (i.e. arrows point one way)
        has_edge(dg, Z, Y) && continue
        inneighbors_Y = inneighbors(dg, Y)
        inneighbors_Z = inneighbors(dg, Z)
        for X in union(inneighbors_Y, inneighbors_Z)
            # X → Y → Z must be a chain (i.e. arrows point one way)
            has_edge(dg, Y, X) && continue
            if has_edge(dg, X, Y) && has_edge(dg, X, Z) && has_edge(dg, Z, X)
                rem_edge!(dg, Z, X)
                verbose && println("⅂  (Rule 2) Removed $Z → $X (rule 2)")
                n_removed += 1
            end
        end
    end
    return n_removed
end

# Rule 3: (avoid creating cycles or new v-structures)
#    X                 X
#  / | \             / | \
# Y  |  W  becomes  Y  |  W
#  ↘ | ↙             ↘ ↓ ↙
#    Z                 Z
function rule3!(alg::PCRobust, dg::SimpleDiGraph; verbose = false)
    n_removed = 0
    candidate_edges = filter(e -> is_bidirectional(dg, e), edges(dg) |> collect)
    used_edges = Vector{SimpleEdge}()
    for edge in candidate_edges
        edge ∉ used_edges || continue
        X, Z = edge.src, edge.dst

        # Two chains X - Y → Z, X - W → Z exist only if Z has at least 2 incoming neighbors,
        # excluding X.
        inZ = setdiff(inneighbors(dg, Z), X)
        length(inZ) ≥ 2 || continue

        for c in combinations(inZ, 2)
            Y, W = c
            if nonadjacent(dg, Y, W)
                rem_edge!(dg, Z, X)
                verbose && println("⅂  (Rule 3) Removed $Z → $X")
                # Once we've directed an edge, don't consider the same edge or
                # the inverted edge again.
                push!(used_edges, SimpleEdge(X, Z))
                push!(used_edges, SimpleEdge(Z, X))
                n_removed += 2
            end
        end
    end
    return n_removed
end

function is_bidirectional(dg::SimpleDiGraph, e::SimpleEdge)
    return has_edge(dg, e.src, e.dst) && has_edge(dg, e.dst, e.src)
end

function nonadjacent(dg::SimpleDiGraph, v1, v2)
    return !has_edge(dg, v1, v2) && !has_edge(dg, v2, v1)
end
