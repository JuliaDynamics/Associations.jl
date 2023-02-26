export OCE

"""
    OCE <: GraphAlgorithm
    OCE(
        utest::IndependenceTest = SurrogateTest(MIShannon(), KSG2(k = 3, w = 3)),
        ctest::C = LocalPermutationTest(CMIShannon(), FPVP(k = 3, w = 3)),
        τmax::T = 5,
        α = 0.05
    )

The optimal causation entropy (OCE) algorithm for causal discovery (Sun et al.,
2015)[^Sun2015].

## Description

The OCE algorithm has three steps to determine the parents of a variable `xᵢ`.
1. Perform pairwise association tests using `utest` and select the variable `xⱼ(-τ)`
    that has the highest significant (i.e. with associated p-value below `α`)
    association with `xᵢ`.
2. Perform conditional independence tests using `ctest`, finding the parent
    `Pₖ` that has the highest association with `xᵢ` given the already selected parents.
    Repeat until no more variables with significant association are found.
3. Backwards elimination of parents `Pₖ` of `xᵢ` for which `xᵢ ⫫ Pₖ | P - {Pₖ}`,
    where `P` is the set of parent nodes found in the previous steps.

## Returns

When used with [`infer_graph`](@ref), it returns a list of parents for each input variable.
In the future, this will return a labelled, directed graph with all the detected
associations.

[^Sun2015]:
    Sun, J., Taylor, D., & Bollt, E. M. (2015). Causal network inference by optimal
    causation entropy. SIAM Journal on Applied Dynamical Systems, 14(1), 73-106.
"""
Base.@kwdef struct OCE{U, C, T} <: GraphAlgorithm
    utest::U = SurrogateTest(MIShannon(), KSG2(k = 3, w = 3))
    ctest::C = LocalPermutationTest(CMIShannon(), FPVP(k = 3, w = 3))
    τmax::T = 5
    α = 0.05
end

function infer_graph(alg::OCE, x)
    parents = select_parents(alg, x)
    return parents
end

"""
    select_parents(alg::OCE, x)

The parent selection step of the [`OCE`](@ref) algorithm, which identifies the
parents of each `xᵢ ∈ x`, assuming that `x` must be integer-indexable, i.e.
`x[i]` yields the `i`-th variable.
"""
function select_parents(alg::OCE, x; verbose = false)

    # Preliminary parents
    τs = Iterators.flatten([-1:-1:-alg.τmax |> collect for xᵢ in x]) |> collect
    js = Iterators.flatten([fill(i, alg.τmax) for i in eachindex(x)]) |> collect
    embeddings = [genembed(xᵢ, -1:-1:-alg.τmax) for xᵢ in x]
    T = typeof(1.0)
    𝒫s = Vector{Vector{T}}(undef, 0)
    for emb in embeddings
        append!(𝒫s, columns(emb))
    end
    # Find the parents of each variable.
    parents = [select_parents(alg, τs, js, 𝒫s, x, k; verbose) for k in eachindex(x)]

    return parents
end

# A simple struct that stores information about selected parents.
Base.@kwdef mutable struct OCESelectedParents{P, PJ, PT}
    i::Int
    parents::P = Vector{Vector{eltype(1.0)}}(undef, 0)
    parents_js::PJ = Vector{Int}(undef, 0)
    parents_τs::PT = Vector{Int}(undef, 0)
end

function Base.show(io::IO, x::OCESelectedParents)
    s = ["x$(x.parents_js[i])($(x.parents_τs[i]))" for i in eachindex(x.parents)]
    all = "Parents for x$(x.i)(0): $(join(s, ", "))"
    show(io, all)
end

function select_parents(alg::OCE, τs, js, 𝒫s, x, i::Int; verbose = false)
    verbose && println("Finding parents for variable x$i(0)")
    idxs_remaining = 1:length(𝒫s) |> collect
    # Account for the fact that the `𝒫ⱼ ∈ 𝒫s` are embedded. This means that some points are
    # lost from the `xᵢ`s.
    xᵢ = @views x[i][alg.τmax+1:end]
    N = length(τs)
    parents = OCESelectedParents(i = i)

    ###################################################################
    # Forward search
    ###################################################################
    # 1. Can we find a significant pairwise association?
    significant_pairwise = select_first_parent!(parents, idxs_remaining, alg, τs, js, 𝒫s, xᵢ; verbose)
    # 2. Continue until there are no more significant conditional pairwise associations
    significant_cond = true
    k = 0
    verbose && println("Conditional tests")
    while significant_cond
        k += 1
        significant_cond = select_conditional_parent!(parents, idxs_remaining, alg, τs, js, 𝒫s, xᵢ; verbose)
    end

    ###################################################################
    # Backward elimination
    ###################################################################
    bw_significant = true
    k = 0
    M = length(parents.parents)
    verbose && println("Backwards elimination")
    k = 1
    while length(parents.parents) >= 1 && k < length(parents.parents)
        verbose && println("\tk=$k, length(parents) = $(length(parents.parents))")
        bw_significant = backwards_eliminate!(parents, alg, xᵢ, k; verbose)
        if bw_significant
            k += 1
        end
    end
    return parents
end

# Pairwise associations
function select_first_parent!(parents, idxs_remaining, alg, τs, js, 𝒫s, xᵢ; verbose = false)
    M = length(𝒫s)

    # Association measure values and the associated p-values
    Is, pvals = zeros(M), zeros(M)
    for (i, Pj) in enumerate(𝒫s)
        test = independence(alg.utest, xᵢ, Pj)
        Is[i] = test.m
        pvals[i] = pvalue(test)
    end

    if all(pvals .>= alg.α)
        verbose && println("\tCouldn't find any significant pairwise associations.")
        return false
    end
    # Select the variable that has the highest significant association with xᵢ.
    # "Significant" means a p-value strictly less than the significance level α.
    Imax = maximum(Is[pvals .< alg.α])
    idx = findfirst(x -> x == Imax, Is)

    if Is[idx] > 0
        verbose && println("\tFound significant pairwise association with: x$(js[idx])($(τs[idx]))")
        push!(parents.parents, 𝒫s[idx])
        push!(parents.parents_js, js[idx])
        push!(parents.parents_τs, τs[idx])
        deleteat!(idxs_remaining, idx)
        return true
    else
        verbose && println("\tCouldn't find any significant pairwise associations.")
        return false
    end
end

function select_conditional_parent!(parents, idxs_remaining, alg, τs, js, 𝒫s, xᵢ; verbose)
    P = Dataset(parents.parents...)
    M = length(𝒫s)
    Is = zeros(M)
    pvals = zeros(M)
    for (i, Pj) in enumerate(𝒫s)
        test = independence(alg.ctest, xᵢ, Pj, P)
        Is[i] = test.m
        pvals[i] = pvalue(test)
    end
    # Select the variable that has the highest significant association with xᵢ.
    # "Significant" means a p-value strictly less than the significance level α.
    if all(pvals .>= alg.α)
        verbose && println("\tCouldn't find any significant pairwise associations.")
        return false
    end
    Imax = maximum(Is[pvals .< alg.α])
    idx = findfirst(x -> x == Imax, Is)

    if Is[idx] > 0
        verbose && println("\tSignificant conditional association with: x$(js[idx])($(τs[idx]))")
        push!(parents.parents, 𝒫s[idxs_remaining[idx]])
        push!(parents.parents_js, js[idxs_remaining[idx]])
        push!(parents.parents_τs, τs[idxs_remaining[idx]])
        deleteat!(idxs_remaining, idx)
        return true
    else
        verbose && println("\tCouldn't find any significant conditional associations.")
        return false
    end
end

function backwards_eliminate!(parents, alg, xᵢ, k; verbose = false)
    M = length(parents.parents)
    P = parents.parents
    Pj = P[k]
    remaining = Dataset(P...)[:, setdiff(1:M, k)]
    test = independence(alg.ctest, xᵢ, Pj, remaining)
    τ, j = parents.parents_τs[k], parents.parents_js[k]
    I = test.m
    # If p-value >= α, then we can't reject the null, i.e. the statistic I is, in
    # the frequentist hypothesis testingworld, indistringuishable from zero.
    if test.pvalue >= alg.α
        verbose && println("\tEliminating k = $k")
        deleteat!(parents.parents, k)
        deleteat!(parents.parents_js, k)
        deleteat!(parents.parents_τs, k)
        return false
    else
        verbose && println("\tpvalue(test)=$(pvalue(test)) > alg.α = $(alg.α)")
        verbose && println("\tNot eliminating anything")
        return true
    end
end
