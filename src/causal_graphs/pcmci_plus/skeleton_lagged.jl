using ShiftedArrays
using Graphs: add_edge!
using MetaGraphs: MetaDiGraph
using MetaGraphs: set_prop!

function skeleton_lagged(alg::PCMCIPlus, X)
    parents = [skeleton_lagged(alg, X, j) for j in eachindex(X)]
    return parents
end

# Find the parents of variable X[i]
function skeleton_lagged(alg::PCMCIPlus, X, j::Int; verbose = false)
    # Currently considered variable.
    xⱼ = X[j]

    # Contemporaneous variables
    Xis = [X[i] for i in eachindex(X)]

    # Potential parents. This is a vector of tuples (j, τ).
    js_τs = lagged_variables(alg, X)

    # Lazily lagged parent vectors.
    lagged_ts = Dict{Tuple{Int, Int}, CircShiftedArray}()
    for (j, τ) in js_τs
        lagged_ts[(j, τ)] = ShiftedArrays.circshift(X[j], τ)
    end

    Imins = fill(Inf, length(js_τs))

    niters = 0
    p = 0

    while length(js_τs) - 1 ≥ p && niters < 100

        niters += 1
        remove_idxs = Int[]
        for (i, (j, τ)) in enumerate(js_τs)
            xᵢτ = lagged_ts[(j, τ)]
            if p == 0
                testresult = independence(alg.pairwise_test, xᵢτ, xⱼ)
            else
                select_vars = [k for k in 1:p if k ≠ i]
                if isempty(select_vars)
                    #println("No vars to pick from => pairwise")
                    testresult = independence(alg.pairwise_test, xᵢτ, xⱼ)
                else
                    #println("There are variables to pick from => conditional")
                    𝒮 = StateSpaceSet((lagged_ts[(j, τ)] for (j, τ) in js_τs[select_vars])...)
                    testresult = independence(alg.pairwise_test, xᵢτ, xⱼ, 𝒮)
                end
            end
            pval = pvalue(testresult)
            I = point_estimate(testresult)
            Imins[i] = min(abs(I), Imins[i])
            if pval < alg.α
                push!(remove_idxs, i)
            end
        end
        # Remove non-significant
        deleteat!(js_τs, remove_idxs)
        deleteat!(Imins, remove_idxs)

        # Sort by measure, largest to smallest
        sort_idxs = sortperm(Imins, rev = true)
        js_τs = js_τs[sort_idxs]
        p += 1
    end

    return js_τs
end

function lagged_variables(alg::PCMCIPlus, X)
    τs = Iterators.flatten([-1:-1:-alg.τmax |> collect for xⱼ in X]) |> collect
    js = Iterators.flatten([fill(i, alg.τmax) for i in eachindex(X)]) |> collect
    return [(j, τ) for (j, τ) in zip(js, τs)]
end
