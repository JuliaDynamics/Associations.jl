using Graphs.SimpleGraphs: SimpleDiGraph
using Graphs: add_edge!

function skeleton_contemporaneous(alg::PCMCIPlus, X, parents)
    all_vars, d = lagged_and_contemp_variables(alg, X)
    n_vertices = length(all_vars)
    dg = SimpleDiGraph(n_vertices)

    for (i, parentsᵢ) in enumerate(parents)
        for j in eachindex(parents)
            if i != j
                add_undirected_edge!(dg, d[(i, 0)], d[(j, 0)])
            end
        end
        for (j, τ) in parentsᵢ
            add_undirected_edge!(dg, d[(i, 0)], d[(j, τ)])
        end
    end
    return dg
end

function lagged_and_contemp_variables(alg::PCMCIPlus, X)
    τs = Iterators.flatten([0:-1:-alg.τmax |> collect for xⱼ in X]) |> collect
    js = Iterators.flatten([fill(i, alg.τmax + 1) for i in eachindex(X)]) |> collect
    js_τs = [(j, τ) for (j, τ) in zip(js, τs)]

    # Also keep track of which node these lagged variables correspond to
    d = Dict{Tuple{Int, Int}, Int}()
    for (i, jτ_tup) in enumerate(zip(js, τs))
        d[jτ_tup] = i
    end

    return js_τs, d
end

function add_undirected_edge!(dg, src_ix, dst_ix)
    add_edge!(dg, src_ix, dst_ix)
    add_edge!(dg, dst_ix, src_ix)
end


function skeleton_contemporaneous(alg::PCMCIPlus, X, j::Int, parents, dg)
    js_τs = contemp_variables(alg, X, j)

    n_links = length(dg) / 2 # directed graph `dg` has two edges per undirected link
    Imins = fill(Inf, n_links)


end


function contemp_variables(alg::PCMCIPlus, X, j::Int)
    τs = [0 for i in eachindex(X) if j != i]
    js = [0 for i in eachindex(X) if j != i]
    js_τs = [(j, τ) for (j, τ) in zip(js, τs)]
    return js_τs
end
