
using Graphs: SimpleGraphs, SimpleGraph, SimpleDiGraph, Edge, AbstractGraph, AbstractEdge
using Graphs: nv, edges, has_edge, inneighbors

"""
    cpdag(alg::PC, skeleton::SimpleDiGraph, separating_sets::Dict{Edge, Vector{Int}}) → dg::SimpleDiGraph

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
function cpdag(alg::PC, skeleton_graph::SimpleDiGraph,
        separating_sets::Dict{SimpleEdge, Vector{Int}}; verbose = false)
    # Orientiation rules are described in a plethora of books and papers in the literature.
    # I found most of them, including the original paper on the PC algorithm, hard to
    # understand, due to either terse and/or ambiguous language. The best
    # non-ambiguous description of the rules I found were in Kalisch & Bühlmann (2008).
    # These are applied here. No effort has been done to make this as efficient as possible.
    verbose && println("Orienting v-structures...")
    # Convert the skeleton to a directed graph.
    dg = orient_vstructures(alg, skeleton_graph, separating_sets; verbose)
    return dg
    edge_was_redirected = trues(4)
    i = 1
    for i = 1:10
        verbose && println("Applying rules sequentially ... (#$i)")
        # Each application of a rule returns the number of edges that were removed/added.
        # If no edge was removed/added in any step, we're done.
        edge_was_redirected[1] = rule1!(alg, dg; verbose)
        edge_was_redirected[2] = rule2!(alg, dg; verbose)
        edge_was_redirected[3] = rule3!(alg, dg; verbose)
        edge_was_redirected[4] = rule4!(alg, dg; verbose)
        if !any(edge_was_redirected)
            verbose && println("No edge was redirected in this iteration. Stopping.")
            break
        end
        i += 1
    end
    verbose && println("Finished orientiation!")

    return dg
end

# Rule 0. Orient v-structures. `X ↔ Y ↔ Z` becomes `X → Y ← Z` if `Y` is not in the
# separating set `S(X, Z)`.
function orient_vstructures(alg::PC, skeleton_graph::SimpleDiGraph, separating_sets::Dict;
        verbose = false)
    # All edges should be
    dg = SimpleDiGraph(skeleton_graph)

    for edge in collect(edges(dg))
        X, Y = edge.src, edge.dst

        Zs = filter(Z -> X != Z && forms_unshielded_triple(dg, X, Y, Z), inneighbors(dg, Y))
        for Z in Zs
            e = SimpleEdge(X, Z)
            re = reverse(e)
            p1 = haskey(separating_sets, e) && Y ∉ separating_sets[e]
            p2 = haskey(separating_sets, re) && Y ∉ separating_sets[re]
            if p1 || p2
                rem_edge!(dg, Y, X)
                rem_edge!(dg, Y, Z)
                verbose && println("  (Rule 0) Oriented $X ↔ $Y ↔ $Z as $X → $Y ← $Z")
            end
        end
    end
    return dg
end

function forms_unshielded_triple(g, i, j, k)
    return is_bidirectional(g, i, j) &&
        is_bidirectional(g, j, k) &&
        adjacent(g, i, j) &&
        adjacent(g, j, k) &&
        !adjacent(g, i, k)
end

function is_undirected(g::AbstractGraph, e::AbstractEdge)
    if has_edge(g, reverse(e))
        return true
    else
        return false
    end
end
is_undirected(g::AbstractGraph, src::Int, targ::Int) = is_undirected(g, SimpleEdge(src, targ))
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
            verbose && println("  (Rule 1) Oriented $X → $Y ↔ $Z as $X → $Y → $Z")
            n_removed += 1
        end
    end
    return n_removed > 0
end

function rule2!(alg::PC, dg::SimpleDiGraph; verbose = false)
    n_removed = 0

    directed_edges = filter(e -> is_directed(dg, e), collect(edges(dg)))
    for edge in directed_edges
        Y, Z = edge.src, edge.dst
        # Find incoming neighbors `Xs` that form a chain `X → Y → Z`
        Xs_Zadjacent = filter(x -> is_directed(dg, x, Y) && adjacent(dg, x, Z), inneighbors(dg, Y))
        for X in Xs_Zadjacent
            @show X, Y, Z
            rem_edge!(dg, Z, X)
            n_removed += 1
            verbose && println("  (Rule 2) Oriented X → Y → Z ↔ X as X → Y → Z ← X")
        end
    end
    return n_removed > 0
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
        Ws = filter(e -> e.dst == Y && nonadjacent(dg, e.src, Y), directed_edges)
        edges_nonadjacent_to_Z = filter(e -> nonadjacent(dg, e.src, Z), alledges)

        for W in Ws
            f = X -> adjacent_and_undirected(dg, X, Y) &&
                adjacent_and_undirected(dg, X, W) &&
                adjacent(X, Z) && is_bidirectional(dg, X, Z)
            Xs = filter(f, edges_nonadjacent_to_Z)
            for X in Xs
                add_edge!(dg, X, Z)
                n_added += 1
                verbose && println("  (Rule 3) Added $X → $Z")
            end
        end
    end
    return n_added > 0
end

# TODO: finish this
# Rule 4: (avoid creating cycles or new v-structures)
#    X                 X
#  ↙ | ↘             ↙ | ↘
# Y  |  W  becomes  Y  |  W
#  ↘ | ↙             ↘ ↓ ↙
#    Z                 Z
function rule4!(alg::PC, dg::SimpleDiGraph; verbose = false)
    n_added = 0

    alledges = collect(edges(dg))
    directed_edges = filter(e -> is_directed(dg, e), alledges)
    for edge in directed_edges
        Y, Z = edge.src, edge.dst
        Ws = filter(e -> e.dst == Y && nonadjacent(dg, e.src, Y), directed_edges)
        edges_nonadjacent_to_Z = filter(e -> nonadjacent(dg, e.src, Z), alledges)

        for W in Ws
            f = X -> adjacent_and_directed(dg, X, Y) &&
                adjacent_and_directed(dg, X, W)
            Xs = filter(f, edges_nonadjacent_to_Z)
            for X in Xs
                add_edge!(dg, X, Z)
                n_added += 1
                verbose && println("  (Rule 4) Added $X → $Z")
            end
        end
    end
    return n_added > 0
end
export nonadjacent, adjacent
function is_bidirectional(dg::SimpleDiGraph, e::SimpleEdge)
    return has_edge(dg, e) && has_edge(dg, reverse(e))
end

function is_bidirectional(dg::SimpleDiGraph, src::Int, dst::Int)
    return is_bidirectional(dg, SimpleEdge(src, dst))
end

# Nonadjacency is defined in Kalisch & Bühlmann as the absence of an edge in either direction
function nonadjacent(dg, v1, v2)
    return !has_edge(dg, v1, v2) && !has_edge(dg, v2, v1)
end

# Adjacency is defined in Kalisch & Bühlmann as the presence of an edge in either direction
function adjacent(dg, v1, v2)
    return has_edge(dg, v1, v2) || has_edge(dg, v2, v1)
end
adjacent(dg, edge) = adjacent(dg, edge.src, edge.dst)

function adjacent_and_undirected(dg, src, dst)
    return adjacent(dg, src, dst) && is_bidirectional(dg, src, dst)
end

function adjacent_and_directed(dg, src, dst)
    return adjacent(dg, src, dst) && has_edge(dg, src, dst)
end
