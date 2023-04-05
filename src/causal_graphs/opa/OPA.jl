using HypothesisTests: OneSampleZTest, pvalue
using ShiftedArrays
using DelayEmbeddings: estimate_delay
using Statistics: median, mean
using Combinatorics: combinations
using Random: default_rng
import Graphs.SimpleGraphs: SimpleDiGraph

export OPA
# TODO: does one even need to make the embeddings completely symmetric? isn't it enough
# to just test *one* time step in the non-causal direction?
"""
    OPA(; τmax = 1, m::Int = 5, α = 0.05,
        measure_pairwise = MIShannon(),
        measure_cond = CMIShannon(),
        est_pairwise = KSG2(k = 5, w = 1),
        est_cond = MesnerShalisi(k = 5, w = 1),
        n_bootstrap = 3000,
        rng = Random.default_rng(),
        eliminate = true,
        stop_predlag::Int = 0)

The naive variant of the optimal predictive asymmetry (PA) algorithm (`OPA`) for inferring
causal graphs. Performs no pre-checks for the forward measure (i.e. full asymemtry
is computed for all possible variables).

## Keyword arguments

- **`stop_predlag::Int`**: During significance testing loops, stop testing if a negative
    predictive asymmetry is found at a prediction lag `η <= stop_predlag`. PA typically
    drops off rapidly after the first few prediction lags, so this can save some
    computation time, at the expense of potentially higher false negative rates. If
    `stop_predlag <= 0`, then the entire test is always performed (no skipping).

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
4. If `eliminate == true`, then perform a backwards elimination step analogous to the
    elimination step in the [`OCE`](@ref) algorithm. If `eliminate == false`, then
    this step is skipped.

## Returns

When used with [`infer_graph`](@ref), a vector of [`OPASelectedParents`](@ref) is
returned, where the `i`-th element of this vector contains the parents for node `i`.
You can convert this vector to a directed graph by using [`SimpleDiGraph`](@ref),
for example

```julia
x = [rand(100) for i = 1:4]
p = infer_graph(OPA(τmax = 3, m = 6), x)
dg = SimpleDiGraph(p)
collect(edges(dg))
```
"""
Base.@kwdef struct OPA{MP, MC, A, K, EP, EC, N, R, EL} <: GraphAlgorithm
    τmax::Int = 1
    m::Int = 6
    α::A = 0.05
    k::K = 1.0 # exponential decay constant
    measure_pairwise::MP = MIShannon()
    measure_cond::MC = CMIShannon()
    est_pairwise::EP = KSG2(k = 20, w = 1)
    est_cond::EC = MesnerShalisi(k = 20, w = 1)
    n_bootstrap::N = 3000
    # TODO: maximum number of conditional discoveries.
    f::Function = mean
    rng::R = Random.default_rng()
    eliminate::EL = false
    stop_predlag::Int = 1
end

function infer_graph(alg::OPA, x; verbose = false)
    parents = select_parents(alg, x; verbose)
    return parents
end

function infer_graph(alg::OPA, x::AbstractDataset; verbose = false)
    parents = select_parents(alg, columns(x); verbose)
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

"""
    OPASelectedParents(i, js, τs)

Parent node indices `js` and lags `τs` that have been found to influence node `i` using
the [`OPA`](@ref) algorithm.
"""
Base.@kwdef mutable struct OPASelectedParents{PJ, PT}
    i::Int
    js::PJ = Vector{Int}(undef, 0)
    τs::PT = Vector{Int}(undef, 0)
end

function SimpleGraphs.SimpleDiGraph(v::Vector{<:CausalityTools.OPASelectedParents})
    N = length(v)
    g = SimpleDiGraph(N)
    for k = 1:N
        parents = v[k]
        for (j, τ) in zip(parents.js, parents.τs)
            if j != k # avoid self-loops
                add_edge!(g, j, k)
            end
        end
    end
    return g
end

function selected(o::OPASelectedParents)
    js, τs = o.js, o.τs
    @assert length(js) == length(τs)
    return join(["x$(js[i])($(τs[i]))" for i in eachindex(js)], ", ")
end

function Base.show(io::IO, x::OPASelectedParents)
    s = ["x$j($τ)" for (j, τ) in zip(x.js, x.τs)]
    if isempty(x.js)
        all = "Evolution of x$(x.i) affected by ∅"
    else
        all = "Evolution of x$(x.i) affected by $(join(s, ", "))"
    end
    show(io, all)
end

# Potential optimization (not sure if it holds analytically, though)
# Optimize first step by perhaps computing raw mutul information. First,
# find the variable that maximizes MI, then run a cascade of signifiacnce
# tests (using asymmetry) until a significant variable is selected.
# Asymmetry is typically highest at small prediction lags, so it might be beneficial
# to look there first.
function select_parents(alg::OPA, x, i::Int; verbose = false)
    (; τmax, m, α, k, measure_pairwise, measure_cond, est_pairwise, est_cond, n_bootstrap,
        f, rng, eliminate, stop_predlag) = alg

    verbose && println("\n=======================================")
    verbose && println("Finding parents for variable x$i ....")
    verbose && println("=======================================")

    # Nodes can be self-causal
    idxs = 1:length(x)
    Pjs = Int[]
    Pτs = Int[]
    for j in idxs
        if j == i
            append!(Pτs, 0:τmax)
            append!(Pjs, repeat([j], τmax + 1))
        else
            append!(Pτs, 0:τmax)
            append!(Pjs, repeat([j], τmax + 1))
        end
    end

    # These are modified in-place during the selection process.
    parents = OPASelectedParents(i = i)

    ###################################################################
    # Forward search
    ###################################################################
    # 1. Can we find a significant pairwise association?
    verbose && println("˧ Querying for pairwise directional association...")

    significant_pairwise = select_first_parent!(alg, parents, x, i, Pjs, Pτs; verbose)
    if significant_pairwise
        verbose && println("˧ Querying new variables conditioned on already selected variables...")
        # 2. Continue until there are no more significant conditional pairwise associations
        significant_cond = true
        k = 0
        while significant_cond && !isempty(Pjs)
            k += 1
            significant_cond = select_conditional_parent!(alg, parents, x, i, Pjs, Pτs; verbose)
        end

        if eliminate
            verbose && println("˧ Backwards elimination...")

            idxs_vars_remaining = 1:length(parents.js) |> collect
            n_initial = length(idxs_vars_remaining)
            q = 0
            while length(idxs_vars_remaining) >= 2 && q < n_initial
                q += 1
                backwards_eliminate!(alg, parents, x, i, q, idxs_vars_remaining; verbose)
            end
        end
    end

    return parents
end

# Strategy
# Find the variable that maximizes MI
# Then perform asymmetry-based significance testing to find the first variable
# that has significant.
function select_first_parent!(alg::OPA, parents, x, i::Int, js, τs; verbose = false)
    (; τmax, m, α, k, measure_pairwise, measure_cond, est_pairwise, est_cond, n_bootstrap,
        f, rng, eliminate, stop_predlag) = alg
    xᵢ = x[i]
    N = length(xᵢ)

    # Generate all candidates
    @assert length(τs) == length(js)

    Δs = [zeros(m) for i in eachindex(js)] # allocating outside is fine, because we overwrite.
    fws = [zeros(m) for i in eachindex(js)] # allocating outside is fine, because we overwrite.

    # Account for circularly shifting the data.
    maxlag = maximum([parents.τs; alg.m])
    idxs = (maxlag + 1):(N - maxlag)
    Ylagged = [lag_for_asymmetry(xᵢ, [0, abs(i)]) for i in 1:m]

    for (ix, (j, τ)) in enumerate(zip(js, τs))
        # Compute asymmetry for the variable with currently highest association with the target.
        X⁻, X⁺ = lag_for_asymmetry(x[j], τ)
        for i in 1:m
            Yⁿ⁻, Yⁿ⁺ = Ylagged[i]
            fw = @views dispatch(measure_pairwise, est_pairwise, X⁻[idxs], Yⁿ⁺[idxs])
            bw = @views dispatch(measure_pairwise, est_pairwise, X⁺[idxs], Yⁿ⁻[idxs])
            if stop_predlag > 0 && i <= stop_predlag && fw - bw < 0
                Δs[ix] .= 0.0
                fws[ix] .= 0.0
                break
            end
            Δs[ix][i] = fw - bw
            fws[ix][i] = fw
        end
        #Δs[ix] = scale_exponentially(Δs[ix]; k)
    end

    # Select the variable that has the highest significant association with xᵢ.
    # "Significant" means a p-value strictly less than the significance level α,
    # and simultanously positive mean forward measure.
    test_results = [bootstrap_right(f, Δ, 0.0; n = n_bootstrap, rng) for Δ in Δs]
    pvals = pvalue.(test_results)
    Is = [mean(Δ) for Δ in Δs]
    Is_fw = [mean(fw) for fw in fws]
    condition = pvals .< alg.α .&& Is_fw .> 0.0 # todo add keyword for latter case

    if !any(condition)
        s = ["x$i(0) ⇍ x($j) | ∅" for j in filter(j -> j != i, js)]
        verbose && println("  $(join(s, "\n  "))")
        return false
    end
    Imax = maximum(Is[condition])
    idx_max = findfirst(x -> x == Imax, Is)

    if !isnothing(idx_max)
        verbose && println("  x$i(0) ← x$(js[idx_max])($(τs[idx_max]))")
        push!(parents.js, js[idx_max])
        push!(parents.τs, τs[idx_max])
        deleteat!(js, idx_max)
        deleteat!(τs, idx_max)
        return true
    end

    # We found no significant pairwise associations
    s = ["x$i(0) ⇍ x($j) | ∅" for j in filter(j -> j != i, js)]
    verbose && println("  $(join(s, "\n  "))")
    return false
end

function select_conditional_parent!(alg::OPA, parents, x, i::Int, js, τs; verbose = false)
    (; τmax, m, α, k, measure_pairwise, measure_cond, est_pairwise, est_cond, n_bootstrap,
        f, rng, eliminate, stop_predlag) = alg

    @assert length(parents.js) == length(parents.τs)
    xᵢ = x[i]
    N = length(xᵢ)

    # Account for circularly shifting the data.
    maxlag = maximum([parents.τs; alg.m])
    idxs = (maxlag + 1):(N - maxlag)
    C⁻, C⁺ = lag_for_asymmetry(x, parents.τs, parents.js)

    Δs = [zeros(m) for i in eachindex(js)] # allocating outside is fine, because we overwrite.
    fws = [zeros(m) for i in eachindex(js)] # allocating outside is fine, because we overwrite.
    Ylagged = [lag_for_asymmetry(xᵢ, [0, abs(i)]) for i = 1:m]
    for (ix, (j, τ)) in enumerate(zip(js, τs))
        # Compute asymmetry for the variable with currently highest association with the target.
        X⁻, X⁺ = lag_for_asymmetry(x[j], τ)
        for i in 1:m
            Yⁿ⁻, Yⁿ⁺ = Ylagged[i]
            fw = @views dispatch(measure_cond, est_cond, X⁻[idxs], Yⁿ⁺[idxs], C⁻[idxs])
            bw = @views dispatch(measure_cond, est_cond, X⁺[idxs], Yⁿ⁻[idxs], C⁺[idxs])
            if stop_predlag > 0 && i <= stop_predlag && fw - bw < 0
                Δs[ix] .= 0.0
                fws[ix] .= 0.0
                break
            end
            Δs[ix][i] = fw - bw
            fws[ix][i] = fw
        end
        #Δs[ix] = scale_exponentially(Δs[ix]; k)
    end

    # Select the variable that has the highest significant association with xᵢ.
    # "Significant" means a p-value strictly less than the significance level α,
    # and simultanously positive mean forward measure.
    test_results = [bootstrap_right(f, Δ, 0.0; n = n_bootstrap, rng) for Δ in Δs]
    pvals = pvalue.(test_results)
    Is = [mean(Δ) for Δ in Δs]
    Is_fw = [mean(fw) for fw in fws]

    condition = pvals .< alg.α .&& Is_fw .> 0.0 # todo add keyword for latter case
    if !any(condition)
        s = ["x$i(0) ⇍ x($j) | ∅" for j in filter(j -> j != i, js)]
        verbose && println("  $(join(s, "\n  "))")
        return false
    end

    Imax = maximum(Is[condition])
    idx_max = findfirst(x -> x == Imax, Is)
    if !isnothing(idx_max)
        verbose && println("  x$i(0) ← x$(js[idx_max])($(τs[idx_max]))")
        push!(parents.js, js[idx_max])
        push!(parents.τs, τs[idx_max])
        deleteat!(js, idx_max)
        deleteat!(τs, idx_max)
        return true
    end

    # No significant links were found.
    s = ["x$i(1) ⇍ x$j($τ) | $(selected(parents)))" for (τ, j) in zip(τs, js)]
    verbose && println("  $(join(s, "\n  "))")
    return false
end

function backwards_eliminate!(alg::OPA, parents, x, i::Int, q::Int, idxs_vars_remaining;
        verbose = false)
    (; τmax, m, α, k, measure_pairwise, measure_cond, est_pairwise, est_cond, n_bootstrap,
        f, rng, eliminate, stop_predlag) = alg

    M = length(parents.js)
    js_remaining = parents.js[setdiff(1:M, q)]
    τs_remaining = parents.τs[setdiff(1:M, q)]
    if isempty(js_remaining) || q > length(parents.js)
        return false
    end
    xᵢ = x[i]
    # Account for circularly shifting the data.
    maxlag = maximum([parents.τs; m])
    idxs = (maxlag + 1):(length(xᵢ) - maxlag)

    # See if x becomes independent of the target, given some subset of the other parents,
    # on conditioning sets of gradually increasing size, starting from 1.
    X⁻, X⁺ = lag_for_asymmetry(x[parents.js[q]], parents.τs[q])

    max_combolength = length(js_remaining)
    Ylagged = [lag_for_asymmetry(xᵢ, [0, abs(i)]) for i = 1:m] # the unlagged xᵢ may now be a parent, so don't include it

    # Go through conditioning sets of increasing size.
    for cl in max_combolength:max_combolength # actually, for now, just do maximum size.
        combs = combinations(1:max_combolength, cl) |> collect # returns combinations of increasing size

        for (k, comb) in enumerate(combs)
            Δs = zeros(m) # allocating outside is fine, because we overwrite.
            fws = zeros(m) # allocating outside is fine, because we overwrite.

            C⁻, C⁺ = lag_for_asymmetry(x, τs_remaining[comb], js_remaining[comb])
            for i in 1:m
                Yⁿ⁻, Yⁿ⁺ = Ylagged[i]
                fw = @views dispatch(measure_cond, est_cond, X⁻[idxs], Yⁿ⁺[idxs], C⁻[idxs])
                bw = @views dispatch(measure_cond, est_cond, X⁺[idxs], Yⁿ⁻[idxs], C⁺[idxs])
                if stop_predlag > 0 && i <= stop_predlag && fw - bw < 0
                    Δs[ix] .= 0.0
                    fws[ix] .= 0.0
                    break
                end
                Δs[i] = fw - bw
                fws[i] = fw
            end
            #Δs = scale_exponentially(Δs; k)
            # If p-value >= α, then we can't reject the null, i.e. the statistic I is
            # indistinguishable from zero, so we claim independence.
            result = bootstrap_right(f, Δs, 0.0; tail = :right)
            if verbose
                if pvalue(result) >= α || mean(fws) <= 0
                    j = parents.js[q]
                    τ = parents.τs[q]
                    r = "Removing x$j($τ) from parent set"
                    s = join(["x$j($τ)" for (j, τ) in zip(js_remaining[comb], τs_remaining[comb])], ", ")
                    verbose && println("  x$i(0) ⫫ x$j($τ) | $s ... $r")
                else
                    j = parents.js[q]
                    τ = parents.τs[q]
                    s = join(["x$j($τ)" for (j, τ) in zip(js_remaining[comb], τs_remaining[comb])], ", ")
                    r = "Keeping x$j($τ) in parent set"
                    verbose && println("  x$i(0) ← x$j($τ) | $s ... $r")
                end
            end
            # A parent became independent of the target conditional on the remaining parents
            if pvalue(result) >= α || mean(fws) <= 0
                deleteat!(parents.js, q)
                deleteat!(parents.τs, q)
                deleteat!(idxs_vars_remaining, q)
                return true
            end
        end
    end
    return false
end
