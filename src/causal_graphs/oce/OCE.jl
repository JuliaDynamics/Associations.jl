export OCE

"""
    OCE <: GraphAlgorithm
    OCE(; utest::IndependenceTest = SurrogateTest(MIShannon(), KSG2(k = 3, w = 3)),
          ctest::C = LocalPermutationTest(CMIShannon(), MesnerShalisi(k = 3, w = 3)),
          τmax::T = 1, α = 0.05)

The optimal causation entropy (OCE) algorithm for causal discovery (Sun et al.,
2015)[^Sun2015].

## Description

The OCE algorithm has three steps to determine the parents of a variable `xᵢ`.
1. Perform pairwise independence tests using `utest` and select the variable `xⱼ(-τ)`
    that has the highest significant (i.e. with associated p-value below `α`)
    association with `xᵢ(0)`. Assign it to the set of selected parents `P`.
2. Perform conditional independence tests using `ctest`, finding the parent
    `Pₖ` that has the highest association with `xᵢ` given the already selected parents,
    and add it to `P`.
    Repeat until no more variables with significant association are found.
3. Backwards elimination of parents `Pₖ` of `xᵢ(0)` for which `xᵢ(0) ⫫ Pₖ | P - {Pₖ}`,
    where `P` is the set of parent nodes found in the previous steps.

`τmax` indicates the maximum lag `τ` between the target variable `xᵢ(0)` and
its potential parents `xⱼ(-τ)`. Sun et al. 2015's method is based on `τmax = 1`.

## Returns

When used with [`infer_graph`](@ref), it returns a vector `p`, where `p[i]` are the
parents for each input variable. In the future, this will return a labelled, directed
graph with all the detected associations.

## Examples

- [Inferring time series graph from a chain of logistic maps](@ref oce_example)

[^Sun2015]:
    Sun, J., Taylor, D., & Bollt, E. M. (2015). Causal network inference by optimal
    causation entropy. SIAM Journal on Applied Dynamical Systems, 14(1), 73-106.
"""
Base.@kwdef struct OCE{U, C, T} <: GraphAlgorithm
    utest::U = SurrogateTest(MIShannon(), KSG2(k = 3, w = 3), nshuffles = 100)
    ctest::C = LocalPermutationTest(CMIShannon(), MesnerShalisi(k = 3, w = 3), nshuffles = 100)
    τmax::T = 1
    α = 0.05
end

function infer_graph(alg::OCE, x; verbose = true)
    parents = select_parents(alg, x; verbose)
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

function selected(o::OCESelectedParents)
    js, τs = o.parents_js, o.parents_τs
    @assert length(js) == length(τs)
    return join(["x$(js[i])($(τs[i]))" for i in eachindex(js)], ", ")
end


function Base.show(io::IO, x::OCESelectedParents)
    s = ["x$(x.parents_js[i])($(x.parents_τs[i]))" for i in eachindex(x.parents)]
    all = "x$(x.i)(0) ← $(join(s, ", "))"
    show(io, all)
end

function select_parents(alg::OCE, τs, js, 𝒫s, x, i::Int; verbose = false)
    verbose && println("\nInferring parents for x$i(0)...")
    # Account for the fact that the `𝒫ⱼ ∈ 𝒫s` are embedded. This means that some points are
    # lost from the `xᵢ`s.
    xᵢ = @views x[i][alg.τmax+1:end]
    N = length(τs)
    parents = OCESelectedParents(i = i)

    ###################################################################
    # Forward search
    ###################################################################
    # 1. Can we find a significant pairwise association?
    verbose && println("˧ Querying pairwise associations...")

    significant_pairwise = select_first_parent!(parents, alg, τs, js, 𝒫s, xᵢ, i; verbose)

    if significant_pairwise
        verbose && println("˧ Querying new variables conditioned on already selected variables...")
        # 2. Continue until there are no more significant conditional pairwise associations
        significant_cond = true
        k = 0
        while significant_cond
            k += 1
            significant_cond = select_conditional_parent!(parents, alg, τs, js, 𝒫s, xᵢ, i; verbose)
        end

        ###################################################################
        # Backward elimination
        ###################################################################
        if !(length(parents.parents) >= 2)
            return parents
        end

        verbose && println("˧ Backwards elimination...")

        eliminate = true
        ks_remaining = Set(1:length(parents.parents))
        while eliminate && length(ks_remaining) >= 2
            for k in ks_remaining
                eliminate = backwards_eliminate!(parents, alg, xᵢ, k; verbose)
                if eliminate
                    filter!(x -> x == k, ks_remaining)
                end
            end
        end
    end
    return parents
end

# Pairwise associations
function select_first_parent!(parents, alg, τs, js, 𝒫s, xᵢ, i; verbose = false)
    M = length(𝒫s)

    if isempty(𝒫s)
        return false
    end

    # Association measure values and the associated p-values
    Is, pvals = zeros(M), zeros(M)
    for (i, Pj) in enumerate(𝒫s)
        test = independence(alg.utest, xᵢ, Pj)
        Is[i] = test.m
        pvals[i] = pvalue(test)
    end

    if all(pvals .>= alg.α)
        s = ["x$i(0) ⫫ x$j(t$τ) | ∅)" for (τ, j) in zip(τs, js)]
        verbose && println("\t$(join(s, "\n\t"))")
        return false
    end
    # Select the variable that has the highest significant association with xᵢ.
    # "Significant" means a p-value strictly less than the significance level α.
    Imax = maximum(Is[pvals .< alg.α])
    idx = findfirst(x -> x == Imax, Is)

    if Is[idx] > 0
        verbose && println("\tx$i(0) !⫫ x$(js[idx])($(τs[idx])) | ∅")
        push!(parents.parents, 𝒫s[idx])
        push!(parents.parents_js, js[idx])
        push!(parents.parents_τs, τs[idx])
        deleteat!(𝒫s, idx)
        deleteat!(js, idx)
        deleteat!(τs, idx)
        return true
    else
        s = ["x$i(0) ⫫ x$j($τ) | ∅)" for (τ, j) in zip(τs, js)]
        verbose && println("\t$(join(s, "\n\t"))")
        return false
    end
end

function select_conditional_parent!(parents, alg, τs, js, 𝒫s, xᵢ, i; verbose)
    if isempty(𝒫s)
        return false
    end

    P = StateSpaceSet(parents.parents...)
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
        s = ["x$i(0) ⫫ x$j($τ) | $(selected(parents))" for (τ, j) in zip(τs, js)]
        verbose && println("\t$(join(s, "\n\t"))")
        return false
    end
    Imax = maximum(Is[pvals .< alg.α])
    idx = findfirst(x -> x == Imax, Is)

    if Is[idx] > 0
        τ = τs[idx]
        j = js[idx]
        verbose && println("\tx$i(0) !⫫ x$j($τ) | $(selected(parents))")
        push!(parents.parents, 𝒫s[idx])
        push!(parents.parents_js, js[idx])
        push!(parents.parents_τs, τs[idx])
        deleteat!(𝒫s, idx)
        deleteat!(τs, idx)
        deleteat!(js, idx)
        return true
    else
        s = ["x$i(1) ⫫ x$j($τ) | $(selected(parents)))" for (τ, j) in zip(τs, js)]
        verbose && println("\t$(join(s, "\n\t"))")
        return false
    end
end

function backwards_eliminate!(parents, alg, xᵢ, k; verbose = false)
    M = length(parents.parents)
    P = parents.parents
    Pj = P[k]
    remaining_idxs = setdiff(1:M, k)
    remaining = StateSpaceSet(P...)[:, remaining_idxs]
    test = independence(alg.ctest, xᵢ, Pj, remaining)

    if verbose
        τ, j = parents.parents_τs[k], parents.parents_js[k] # Variable currently considered
        τs = parents.parents_τs
        js = parents.parents_js
        src_var = "x$j($τ)"
        targ_var = "x$(js[k])($(τs[k]))"
        cond_var = join(["x$(js[i])($(τs[i]))" for i in remaining_idxs], ", ")

        if test.pvalue >= alg.α
            outcome_msg = "Removing x$(j)($τ) from parent set"
            println("\t$src_var ⫫ $targ_var | $cond_var → $outcome_msg")
        else
            outcome_msg = "Keeping x$(j)($τ) in parent set"
            println("\t$src_var !⫫ $targ_var | $cond_var → $outcome_msg")
        end
    end

    # If p-value >= α, then we can't reject the null, i.e. the statistic I is
    # indistinguishable from zero, so we claim independence and remove the variable.
    if test.pvalue >= alg.α
        deleteat!(parents.parents, k)
        deleteat!(parents.parents_js, k)
        deleteat!(parents.parents_τs, k)
        return true
    else
        return false
    end
end
