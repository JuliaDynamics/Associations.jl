
using Graphs: SimpleGraphs, SimpleGraph, SimpleDiGraph, Edge, AbstractGraph, AbstractEdge
using Graphs: nv, edges, has_edge, inneighbors

"""
    cpdag(alg::PC, skeleton::SimpleGraph, separating_sets::Dict{Edge, Vector{Int}}) → dg::SimpleDiGraph

Orient edges in the `skeleton` graph using the given `separating_sets` using
algorithm 2 in Kalisch & Bühlmann (2008[^Kalisch2008]) and return the directed graph `cpdag`.

## Description

First, a directed graph `dg` is constructed from `skeleton`, such that every
undirected edge `X - Y` in `skeleton` is replaced by the bidirectional edge `X ↔ Y` in `dg`.
In practices, for each `X ↔ Y`, we construct two directional edges `X → Y` and `Y → X`.

Orientiation rules 0-3 are then applied to `dg`. We use the rules as stated in
Colombo & Maathuis, 2014.
- Rule 1 (prevents new v-structures): `X → Y ↔ Z` becomes `X → Y → Z` if `X` and `Z`
    are not adjacent.
- Rule 2 (avoids cycles): `X → Y → Z ↔ X` becomes X → Y → Z ← X`
- Rule 3: To avoid creating cycles or new v-structures,
       X                 X
     | | |             | | |
    Y  |  W  becomes  Y  |  W
     ↘ | ↙             ↘ ↓ ↙
       Z                 Z


- Rule 0 (orients v-structures): `X ↔ Y ↔ Z` becomes `X → Y ← Z` if `Y` is not in the
    separating set `S(X, Z)`.


[^Kalisch2008]:
    Kalisch, M., & Bühlmann, P. (2008). Robustification of the PC-algorithm for directed
    acyclic graphs. Journal of Computational and Graphical Statistics, 17(4), 773-789.
"""
function cpdag(alg::PC, skeleton::SimpleDiGraph,
        separating_sets::Dict{SimpleEdge, Vector{Int}}; verbose = false)
    # Convert the skeleton to a directed graph where (Xᵢ - Xⱼ) becomes (Xᵢ → Xⱼ, Xⱼ → Xᵢ).
    dg = deepcopy(skeleton)

    # Orientiation rules are described in a plethora of books and papers in the literature.
    # I found most of them, including the original paper on the PC algorithm, hard to
    # understand, due to either terse and/or ambiguous language. The best
    # non-ambiguous description of the rules I found were in Kalisch & Bühlmann (2008).
    # These are applied here. No effort has been done to make this as efficient as possible.
    verbose && println("Transforming skeleton graph into cpdag...")
    verbose && println("⅂ Level 0:")
    rule0!(alg, dg, separating_sets; verbose)

    edge_was_redirected = ones(3)
    i = 1
    while any(edge_was_redirected .> 0) && i < alg.maxiters_orient
        verbose && println("⅂ Level 1:")
        # Each application of a rule returns the number of edges that were removed/added.
        # If no edge was removed/added in any step, we're done.
        edge_was_redirected[1] = rule1!(alg, dg; verbose)
        edge_was_redirected[2] = rule2!(alg, dg; verbose)
        edge_was_redirected[3] = rule3!(alg, dg; verbose)
        i += 1
    end
    return dg
end

# Rule 0. Orient v-structures. `X ↔ Y ↔ Z` becomes `X → Y ← Z` if `Y` is not in the
# separating set `S(X, Z)`.
function rule0!(alg::PC, dg::SimpleDiGraph, separating_sets::Dict;
        verbose = false)
    n_removed = 0
    for edge in edges(dg)
        verbose && println("Considering $edge")
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
    return n_removed
end

function is_undirected(g::AbstractGraph, e::AbstractEdge)
    if has_edge(g, reverse(e))
        return true
    else
        return false
    end
end

function is_undirected(g::AbstractGraph, src::Int, targ::Int)
    e = SimpleEdge(src, targ)
    if has_edge(g, reverse(e))
        return true
    else
        return false
    end
end

is_directed(args...) = !is_undirected(args...)

# TODO: is the fact that we're directly modifying the graph a problem? Should we make a
# copy?
function rule1!(alg::PC, dg::SimpleDiGraph; verbose = false)
    n_removed = 0

    directed_edges = filter(e -> is_directed(dg, e), collect(edges(dg)))
    for edge in directed_edges
        X::Int = edge.src
        Y::Int = edge.dst

        # Find the node numbers `Zs` of all incoming undirected edges to Y, except `X`,
        # which is already part of `edge`.
        undirected_Zs_nonadjacent_to_X = setdiff(inneighbors(dg, Y), X)
        filter!(Z -> is_undirected(dg, Z, Y) && nonadjacent(dg, Z, X), undirected_Zs_nonadjacent_to_X)

        for Z in undirected_Zs_nonadjacent_to_X
            rem_edge!(dg, Z, Y) # Direct the undirected edge
            verbose && println("⅂  (Rule 1) Oriented $X → $Y ↔ $Z as $X → $Y → $Z")
            n_removed += 1
            return n_removed
        end
    end
    return n_removed
end

function rule2!(alg::PC, dg::SimpleDiGraph; verbose = false)
    n_removed = 0

    directed_edges = filter(e -> is_directed(dg, e), collect(edges(dg)))
    for edge in directed_edges
        Y, Z = edge.src, edge.dst

        # Find incoming neighbors `Xs` that form a chain `X → Y → Z`
        Xs_Zadjacent = filter(x -> is_directed(dg, x, Y) && adjacent(dg, x, Z), inneighbors(dg, Y))

        for X in Xs_Zadjacent
            # X → Y → Z must be a chain (i.e. arrows point one way)
            rem_edge!(dg, Z, X)
            n_removed += 1
            verbose && println("⅂  (Rule 2) Oriented X → Y → Z ↔ X as X → Y → Z ← X")
            return n_removed
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
function rule3!(alg::PC, dg::SimpleDiGraph; verbose = false)
    n_added = 0

    alledges = collect(edges(dg))
    directed_edges = filter(e -> is_directed(dg, e), alledges)
    for edge in directed_edges
        Y, Z = edge.src, edge.dst
        Ws = filter(e -> e.dst == Y, directed_edges)
        edges_nonadjacent_to_Z = filter(e -> nonadjacent(dg, e.src, Z), alledges)

        for W in Ws
            f = X -> adjacent_and_undirected(dg, X, Y) && adjacent_and_undirected(dg, X, W)
            Xs = filter(f, edges_nonadjacent_to_Z)
            for X in Xs
                add_edge!(dg, X, Z)
                n_added += 1
                return n_added
                verbose && println("⅂  (Rule 3) Added $X → $Z")
            end
        end
    end
    return n_added
end

function is_bidirectional(dg::SimpleDiGraph, e::SimpleEdge)
    return has_edge(dg, e.src, e.dst) && has_edge(dg, e.dst, e.src)
end

function is_bidirectional(dg::SimpleDiGraph, src::Int, dst::Int)
    return has_edge(dg, src, dst) && has_edge(dg, dst, src)
end

# Nonadjacency is defined in Kalisch & Bühlmann as the absence of an edge in either direction
function nonadjacent(dg::SimpleDiGraph, v1, v2)
    return !has_edge(dg, v1, v2) && !has_edge(dg, v2, v1)
end

# Adjacency is defined in Kalisch & Bühlmann as the presence of an edge in either direction
function adjacent(dg::SimpleDiGraph, v1, v2)
    return has_edge(dg, v1, v2) || has_edge(dg, v2, v1)
end

function adjacent(dg::SimpleDiGraph, edge)
    return has_edge(dg, edge) || has_edge(dg, reverse(edge))
end

function adjacent_and_undirected(dg, src, dst)
    return adjacent(dg, src, dst) && is_bidirectional(dg, src, dst)
end
