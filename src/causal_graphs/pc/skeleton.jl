using Graphs: SimpleGraph, SimpleDiGraph, SimpleEdge
using Graphs: nv, complete_digraph, rem_edge!
using Graphs.SimpleGraphs: all_neighbors
using Combinatorics: powerset, combinations

#export skeleton

"""
    skeleton(algorithm::PC, x) → (g, s)

Infer the skeleton graph for the variables `x` using the provided `algorithm`.
`x` must be some iterable where the `i`-th variable can be accessed as `x[i]`.

Return a tuple of the undirected skeleton graph `g::SimpleGraph`, and the
separating sets `s::Dict{SimpleEdge, Vector{Int}})`.
"""
function skeleton(alg::PC, x; verbose = false)
    N = length(x)
    if alg.maxdepth == Inf
        max_degree = N - 2
    else
        max_degree = alg.maxdepth
    end
    graph = complete_digraph(N)
    separating_set = Dict{SimpleEdge, Vector{Int}}()

    skeleton_unconditional!(alg, graph, x; verbose) # only considers pairs
    for 𝓁 = 1:max_degree
        skeleton_conditional!(alg, graph, separating_set, x, 𝓁; verbose)
    end
    return graph, separating_set
end

"""
    skeleton_unconditional!(alg::PC, g::SimpleGraph, x)

Perform pairwise independence tests between all vertices in the `graph`, where
each vertex is represented by the data `x[i]`, and remove any edges in `graph`
where adjacent vertices are found to be independent according to the given independence
`test`. The null hypothesis of independence is rejected
whenever the p-value is below `α`.

This is essentially algorithm 3.2 in Colombo & Maathuis (2014), but only
considering the case `𝓁 = 0`.

Modifies `graph` in-place.

If `alg.pairwise_test` is a directed test, then edges are considered one-by-one.
If `alg.pairwise_test` is not a directed test, then edges (`X → Y, `Y → X`)
are considered simultaneously.

[^Colombo2014]:
    Colombo, D., & Maathuis, M. H. (2014). Order-independent constraint-based causal
    structure learning. J. Mach. Learn. Res., 15(1), 3741-3782.
"""
function skeleton_unconditional!(alg::PC, graph::SimpleDiGraph, x; verbose = false)
    N = length(x)
    pairs = (Tuple(pair) for pair in combinations(1:N, 2))
    for pair in pairs
        s, t = pair
        # If pval > α, then, based on the given the data, we can't reject the hypothesis
        # that `x[s] ⫫ x[t]`. Therefore, we assume that they *are* independent.
        indep_test = independence(alg.pairwise_test, x[s], x[t])
        pval = pvalue(indep_test)
        if pval > alg.α
            edge1 = SimpleEdge(s, t)
            edge2 = SimpleEdge(t, s)
            verbose && println("Skeleton, pairwise: Removing $edge1 and $edge2 (p = $pval)")
            rem_edge!(graph, edge1)
            rem_edge!(graph, edge2)
        end
    end
    return graph
end

"""
    skeleton_conditional!(alg::PC, graph, x, conditional_test::IndependenceTest)

Thin the skeleton `graph`, where each vertex is represented by the data `x[i]`,
by using `conditional_test`. Whenever `x[i] ⫫ x[j] | x[S]` for
some set of variables not including `i` and `j`, the edge between `i` and `j` is
removed, and `S` is stored in the separating set for `i` and `j`.
This is essentially algorithm 3.2 in Colombo & Maathuis (2014), for
the cases `𝓁 >= 1`.

Modifies `graph` in-place.

[^Colombo2014]:
    Colombo, D., & Maathuis, M. H. (2014). Order-independent constraint-based causal
    structure learning. J. Mach. Learn. Res., 15(1), 3741-3782.
"""
function skeleton_conditional!(alg::PC, graph, separating_set, x, 𝓁::Int;
        verbose = false)

    N = length(x)
    # `a[i]` := adjacent vertices to vertex `i`
    a = [all_neighbors(graph, i) for i in 1:nv(graph)]
    ctr = 0
    for (i, aᵢ) in enumerate(a)
        for j in aᵢ
            # The powerset of remaining variables (not including variable i or variable j),
            # limited to subsets of cardinality 𝓁 <= C <= 𝓁 + 1.
            𝐒 = powerset(setdiff(aᵢ, j), 𝓁, 𝓁 + 1) |> collect

            # Perform independence tests and remove the edge between i and j if
            # they are found to be independent.
            ctr += conditionaltest_and_remove_edge!(alg, x, 𝐒, 𝓁, i, j, graph, separating_set; verbose)
        end
    end
    return graph, ctr
end

function conditionaltest_and_remove_edge!(alg::PC, x, 𝐒, 𝓁, i, j, graph, separating_set;
        verbose = false)
    # If there's at least one available variable to condition on.
    ctr = 0

    if length(𝐒) >= 𝓁
        src, trg = @views x[i], x[j]
        # For each subset of variables (not including i and j), perform a conditional
        # independence test `i ⫫ j | Sₖ`. If this holds for any variable(s) `Sₖ`,
        # then the variables are taken as independent, and `Sₖ` is assigned to the
        # separating set for the edges i - j (or, equivalently, j - i).

        # Only pick conditional sets with valid sizes
        if isnothing(alg.maxdepth)
            conditional_sets = 𝐒
        else
            conditional_sets = filter(s -> length(s) <= alg.maxdepth, 𝐒)
        end
        for Sₖ in conditional_sets
            Ŝ = @views StateSpaceSet(x[Sₖ]...)
            # If pval > α, then, based on the given the data, we can't reject the hypothesis
            # that `src ⫫ trg | Ŝ`. Therefore, we assume that they *are* independent
            # given Ŝ.
            res = independence(alg.conditional_test, src, trg, Ŝ)
            pval = pvalue(res)

            if pval > alg.α
                edge1 = SimpleEdge(i, j)
                edge2 = SimpleEdge(j, i)

                verbose && println("Skeleton (conditional, level 𝓁=$𝓁): Removing $edge1 and $edge2 (p = $pval)")
                rem_edge!(graph, edge1)
                rem_edge!(graph, edge2)

                separating_set[edge1] = Sₖ
                separating_set[edge2] = Sₖ
                ctr +=1
                break
            end
        end
    end
    return ctr
end
