using HypothesisTests: SignedRankTest, pvalue
using ShiftedArrays
using DelayEmbeddings: estimate_delay
using Statistics: median

export OPA

"""
    OPA(τmax = 1, m::Int = 5, m = 5, α = 0.05, est = FPVP(), k = 0.5)

The optimal predictive asymmetry (OPA) for inferring causal graphs (Haaga et al., in prep).

## Description

Assumes two inputs: a timeseries ``X^{(i)}``, and a set of candidate timeseries
``X^{(k)}_{k=1}^r`` that potentially influences ``X^{(i)}``.

Under the assumption that the random variable ``X^{(i)}`` is driven by its own past,
the OPA algorithm attempts to quantify if any other input variables ``X^{(j)}``
influences ``X^{(i)}`` by:

1. Computing the pairwise predictive asymmetry distribution with with respect to ``X^{(i)}``
    for all ``X^{(j)} \\in X^{(k)}_{k=1}^r``. Select the variable ``X^{(j)}`` that has
    the highest significant asymmetry  (i.e. the associated p-value is
    below `α`). If an ``X^{(j)}`` is found, it is added to the set of already selected
    parents ``\\mathcal{P}^{(i)}`` for ``X^{(j)}``.
2. Perform conditional OPA, finding the parent
    `Pₖ` that has the highest association with ``X^{(i)}`` given ``X^{(i)}``'s past
    and the already selected parents.
3. Repeat until no more variables with significant association are found.
"""
Base.@kwdef struct OPA{A, K, E}
    τmax::Int = 0
    m::Int = 5
    α::A = 0.05
    k::K = 0.0 # exponential decay constant
    est::E = FPVP(k = 3, w = 3)
end

function infer_graph(alg::OPA, x; verbose = false)
    parents = select_parents(alg, x; verbose)
    return parents
end

"""
    select_parents(alg::OPA, x)

The parent selection step of the [`OPA`](@ref) algorithm, which identifies the
parents of each `xᵢ ∈ x`, assuming that `x` must be integer-indexable, i.e.
`x[i]` yields the `i`-th variable.
"""
function select_parents(alg::OPA, x; verbose = false)
    parents = [select_parents(alg, x, k; verbose) for k in eachindex(x)]
    return parents
end


Base.@kwdef mutable struct OPASelectedParents{P, PJ, PT}
    i::Int
    parents::P = Vector{Vector{eltype(1.0)}}(undef, 0)
    js::PJ = Vector{Int}(undef, 0)
    τs::PT = Vector{Int}(undef, 0)
end

function Base.show(io::IO, x::OPASelectedParents)
    s = ["x$(x.js[i])($(x.τs[i]))" for i in eachindex(x.js)]
    if isempty(x.js)
        all = "Parents for x$(x.i): ∅"
    else
        all = "Parents for x$(x.i): $(join(s, ", "))"
    end
    show(io, all)
end

function select_parents(alg::OPA, x, i::Int; verbose = false)
    verbose && println("Finding parents for variable x$i")

    # These are modified in-place during the selection process.
    js = setdiff(1:length(x), i)
    parents = OPASelectedParents(i = i)

    ###################################################################
    # Forward search
    ###################################################################
    # 1. Can we find a significant pairwise association?
    significant_pairwise = select_first_parent!(alg, parents, x, i, js; verbose)
    significant_conditional = true
    if significant_conditional && significant_pairwise && length(js) >= 1
        # TODO: While loop. This only runs once.
        significant_conditional = select_conditional_parent!(alg, parents, x, i, js; verbose)
    end

    # # TODO: backwards elimination step
    #print("Backwards elimination")

    return parents
end

function scale_exponentially(x, y; k = 1)
    ds = [(η-1) for η in x]
    return y .* exp.(.-ds .* k)
end

# TODO: first test mutual information, then
function select_first_parent!(alg::OPA, parents, x, i::Int, js; verbose = false)
    est = alg.est
    S⁰s = @views [x[j] for j in js]  # nonshifted
    T⁰ = @views x[i] # nonshifted
    Tⁿ⁺s = [ShiftedArrays.circshift(T⁰, -η) for η in 1:alg.m]
    Tⁿ⁻s = [ShiftedArrays.circshift(T⁰, η) for η in 1:alg.m]
    ΔA = [zeros(alg.m) for j in js]

    #τT = estimate_delay(T⁰, "mi_min")
    for (l, j) in enumerate(js)
        ix = findfirst(idx -> idx == j, js)
        S⁰ = S⁰s[ix]
        #τS = estimate_delay(S⁰, "mi_min")

        S⁺ = ShiftedArrays.circshift(S⁰, -1) # TODO: flexible lags here
        S⁻ = ShiftedArrays.circshift(S⁰, 1) # TODO: flexible lags here.
        for (k, η) in enumerate(1:alg.m)
            @views Tⁿ⁺ = Tⁿ⁺s[k]
            @views Tⁿ⁻ = Tⁿ⁻s[k]
            T⁺ = ShiftedArrays.circshift(T⁰, -1) # TODO: flexible lags here
            T⁻ = ShiftedArrays.circshift(T⁰, 1)  # TODO: flexible lags here

            F⁻ = Dataset(Tⁿ⁻)#Dataset(Tⁿ⁻, T⁰, Sⁿ⁻)#, Sⁿ⁻)
            F⁺ = Dataset(Tⁿ⁺)#Dataset(Tⁿ⁺, T⁰, Sⁿ⁺)#, Sⁿ⁺)
            CP⁺ = Dataset(T⁰, T⁺, S⁺)#Dataset(Tⁿ⁺) # chosen parents
            CP⁻ = Dataset(T⁰, T⁻, S⁻)#Dataset(Tⁿ⁻) # chosen parents
            PP⁺ = Dataset(S⁰) # potential parents
            PP⁻ = Dataset(S⁰) # potential parents

            # TODO: there is bias now due to circular shifting. Remove start/end points,
            # preferably using a non-allocating view.
            # A(S → T) = A(S⁰ → T⁺ | T⁰)
            fw = condmutualinfo(CMIShannon(), est, F⁺, PP⁻, CP⁻)
            bw = condmutualinfo(CMIShannon(), est, F⁻, PP⁺, CP⁺)
            #@views fw = condmutualinfo(CMIShannon(), est, S⁰, T⁰, Tⁿ⁺)#  F⁺, PP, CP⁻)
            #@views bw = condmutualinfo(CMIShannon(), est, S⁰, T⁰, Tⁿ⁻)#F⁻, PP, CP⁺)
            ΔA[l][k] = fw - bw
        end
    end
    #@show ΔA[1]
    scaled_ΔA = [scale_exponentially(1:alg.m, a, k = alg.k) for a in ΔA]
    Is = median.(scaled_ΔA)  #mean.(ΔA) #mean.(scaled_ΔA)
    tests = SignedRankTest.(scaled_ΔA) #OneSampleTTest.(ΔA) #OneSampleTTest.(scaled_ΔA)
    pvals = pvalue.(tests, tail = :right)

    if !any(pvals .< alg.α) # No significant links  found.
        s = ["x$i ⫫ x($j) | ∅" for j in filter(j -> j != i, js)]
        verbose && println("\t$(join(s, "\n\t"))")
        return false
    end
    Imax = maximum(Is[pvals .< alg.α])
    idx = findfirst(x -> x == Imax, Is)
    if !isnothing(idx)
        verbose && println("\tx$i ← x$(js[idx])")
        push!(parents.parents, x[idx])
        push!(parents.js, js[idx])
        push!(parents.τs, 0)
        deleteat!(js, idx)
        return true
    end
end

# i is the target variable indexable
function select_conditional_parent!(alg::OPA, parents, x, i::Int, potential_js; verbose = false)
    est = alg.est
    T⁰ = @views x[i] # nonshifted
    Tⁿ⁺s = [ShiftedArrays.circshift(T⁰, -η) for η in 1:alg.m]
    Tⁿ⁻s = [ShiftedArrays.circshift(T⁰, η) for η in 1:alg.m]
    ΔA = [zeros(alg.m) for j in potential_js]

    already_selected = Dataset(parents.parents...)
    already_selected⁻ = Dataset([ShiftedArrays.circshift(x, 1) for x in parents.parents]...)
    already_selected⁺ = Dataset([ShiftedArrays.circshift(x, -1) for x in parents.parents]...)
    #τT = estimate_delay(T⁰, "mi_min")
    T⁺ = ShiftedArrays.circshift(T⁰, -1)# TODO: flexible
    T⁻ = ShiftedArrays.circshift(T⁰, 1) #TODO: flexible

    for (l, j) in enumerate(potential_js)
        S⁰ = @views x[j]
        #τS = estimate_delay(S⁰, "mi_min")
        S⁺ = ShiftedArrays.circshift(S⁰, -1) # TODO: flexible
        S⁻ = ShiftedArrays.circshift(S⁰, 1) # TODO:flexible
        for (k, η) in enumerate(1:alg.m)
            @views Tⁿ⁺ = Tⁿ⁺s[k]
            @views Tⁿ⁻ = Tⁿ⁻s[k]
            F⁻ = Dataset(Tⁿ⁻, T⁰, S⁻) # target past
            F⁺ = Dataset(Tⁿ⁺, T⁰, S⁺) # target future
            CP⁺ = Dataset(T⁺, already_selected⁺, already_selected) # chosen parents future
            CP⁻ = Dataset(T⁻, already_selected⁻, already_selected) # chosen parents past
            PP = Dataset(S⁰) # potential parents

            # TODO: there is bias now due to circular shifting. Remove start/end points,
            # preferably using a non-allocating view.
            # A(S → T) = A(S⁰ → T⁺ | T⁰)
            fw = condmutualinfo(CMIShannon(), est, F⁺, PP, CP⁻)
            bw = condmutualinfo(CMIShannon(), est, F⁻, PP, CP⁺)

            ΔA[l][k] = fw - bw
        end
    end

    scaled_ΔA = [scale_exponentially.(1:alg.m, a, k = alg.k) for a in ΔA]
    Is = median.(scaled_ΔA)  #mean.(ΔA) #mean.(scaled_ΔA)
    tests = SignedRankTest.(scaled_ΔA) #OneSampleTTest.(ΔA) #OneSampleTTest.(scaled_ΔA)
    pvals = pvalue.(tests, tail = :right)

    if !any(pvals .< alg.α) # No significant links found.
        s = ["x$i ⫫ x($j) | x($(parents.js)))" for j in potential_js]
        verbose && println("\t$(join(s, "\n\t"))")
        return false
    end
    Imax = maximum(Is[pvals .< alg.α])
    idx = findfirst(x -> x == Imax, Is)
    if !isnothing(idx)
        verbose && println("\tx$i ← x($(potential_js[idx])) given x($(parents.js)))")
        push!(parents.parents, x[potential_js[idx]])
        push!(parents.js, potential_js[idx])
        push!(parents.τs, 0)
        return true
    end
end
