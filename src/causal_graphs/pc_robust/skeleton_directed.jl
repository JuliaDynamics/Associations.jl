using Graphs: SimpleDiGraph, SimpleEdge
using Graphs: add_edge!, outneighbors
using Combinatorics: permutations

# TODO: this is only designed for non-directed measures. Use type system?
function infer_graph(algorithm::PCRobust{<:IndependenceTest{<:DirectedAssociationMeasure}}, x; verbose = false)
    skeleton_graph, separating_sets = skeleton(algorithm, x; verbose)
    #return skeleton_graph
    #return skeleton_graph = skeleton(algorithm, x; verbose)
    directed_graph = cpdag(algorithm, skeleton_graph, separating_sets; verbose)
    #@show edges(directed_graph) |> collect
    return directed_graph
end


function skeleton(alg::PCRobust{<:IndependenceTest{<:DirectedAssociationMeasure}},
        x; verbose = false)

    N = length(x)
    if alg.maxdepth isa Nothing
        max_degree = N - 2
    else
        max_degree = alg.maxdepth
    end
    graph = SimpleDiGraph(N)
    separating_set = Dict{SimpleEdge, Vector{Int}}()

    skeleton_unconditional!(alg, graph, x; verbose) # only considers pairs
    #@show edges(graph) |> collect
    for 𝓁 = 1:max_degree
        skeleton_conditional!(alg, graph, separating_set, x, 𝓁; verbose)
    end
    return graph, separating_set
end


function skeleton_unconditional!(
        alg::PCRobust{<:IndependenceTest{<:DirectedAssociationMeasure}},
        graph::SimpleDiGraph,
        x;
        verbose = false)
    N = length(x)
    @show "directed unconditional"
    pairs = (Tuple(pair) for pair in permutations(1:N, 2))
    for pair in pairs
        s, t = pair
        # If pval > α, then, based on the given the data, we can't reject the hypothesis
        # that `x[s] ⫫ x[t]`. Therefore, we assume that they *are* independent.
        pval = @views pvalue(independence(alg.unconditional_test, x[s], x[t]))
        if pval < alg.α
            @show s, t, pval, alg.α
            edge = SimpleEdge(s, t)
            verbose && println("Skeleton, pairwise: Adding $edge (p = $pval)")
            add_edge!(graph, edge)
        end
    end

    return graph
end

function skeleton_conditional!(
        alg::PCRobust{<:IndependenceTest{<:DirectedAssociationMeasure}},
        graph,
        separating_set, x, 𝓁::Int;
        verbose = false)
    @show "directed conditional"

    N = length(x)
    # `a[i]` := incoming vertices to vertex `i`, i.e. only vertices that have
    # been found to have a direct connection to  vertex `i`
    a = [outneighbors(graph, i) for i in 1:nv(graph)]
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

function conditionaltest_and_remove_edge!(
        alg::PCRobust{<:IndependenceTest{<:DirectedAssociationMeasure}},
        x, 𝐒, 𝓁, i, j, graph, separating_set;
        verbose = false)
    # If there's at least one available variable to condition on.
    ctr = 0
    if length(𝐒) >= 𝓁
        s, t = @views x[i], x[j]
        # For each subset of variables (not including i and j), perform a conditional
        # independence test `i ⫫ j | Sₖ`. If this holds for any variable(s) `Sₖ`,
        # then the variables are taken as independent, and `Sₖ` is assigned to the
        # separating set for the edges i - j (or, equivalently, j - i).

        # Only pick conditional sets with valid sizes
        if isnothing(alg.maxdepth)
            conditional_sets = 𝐒
        else
            conditional_sets = 𝐒[length.(𝐒) .<= alg.maxdepth]
        end
        for Sₖ in conditional_sets
            Ŝ = @views Dataset(x[Sₖ]...)
            # If pval > α, then, based on the given the data, we can't reject the hypothesis
            # that `src ⫫ trg | Ŝ`. Therefore, we assume that they *are* independent
            # given Ŝ.

            res = independence(alg.conditional_test, s, t, Ŝ)
            pval = pvalue(res)

            if pval >= alg.α
                @show i, j, Sₖ, pval, alg.α

                edge = SimpleEdge(i, j)
                verbose && println("Skeleton (conditional, level 𝓁=$𝓁): Removing $edge (p = $pval)")
                rem_edge!(graph, edge)
                separating_set[edge] = Sₖ
                ctr +=1
                break
            end
        end
    end
    return ctr
end
