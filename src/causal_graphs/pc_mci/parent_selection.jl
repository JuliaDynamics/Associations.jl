export select_parents

"""
    select_parents(alg::PCMCI, x, i::Int)

The parent selection step of the [`PCMCI`](@ref) algorithm, which identifies the
parents of each `xᵢ ∈ x`, assuming that `x` must be integer-indexable, i.e.
`x[i]` yields the `i`-th variable.

This is algorithm S1 in Runge et al. (2018).
"""
function select_parents(alg::PCMCI, x)
    # Preliminary parents
    τs = Iterators.flatten([-1:-1:-alg.τmax |> collect for xᵢ in x]) |> collect
    embeddings = [columns(genembed(xᵢ, τs[i])) for (i, xᵢ) in enumerate(x)]
    𝒫s = Vector(undef, 0)
    for emb in embeddings
        append!(𝒫s, columns(emb))
    end
    return select_parents(alg, 𝒫s, τs, x, 1)
end

"""
    select_parents!(alg::PCMCI, 𝒫s::Vector{AbstractStateSpaceSet}, x, j::Int)

Select parents for variable `x[j]`, given the embeddings `𝒫s`, where `𝒫s[j]` is the
`alg.τmax`-dimensional embedding for `x[j]` where columns correspond to
embeddings lags `1:alg.τmax`.
"""
function select_parents(alg::PCMCI, 𝒫s, τs, x, j)
    q = alg.
    Is = [fill(Inf, alg.τmax) for xᵢ in x]
    pvals = [fill(1.0, alg.τmax) for xᵢ in x]

    js = [repeat([i], alg.τmax) for xᵢ  in x]

     ## Preliminary parents.
     𝒫 = Iterators.flatten(𝒫s)

    # Keeps track of which variables are remaining.
    τsⱼ = Iterators.flatten(copy(τs))


    return 𝒫

    for D = 1:alg.pmax # loop over conditioning set dimension.
        # Break if there are no more variables to select.
        n_unselected = sum(length.(τsⱼ))
        if n_unselected < cond_dim
            break
        end


        for (i, Pᵢ) in enumerate(𝒫) # Line 9
            q = 0
            # For all subsets 𝒮 ⊆ (𝒫 excluding 𝒫[i]) s.t. card(𝒮) = D
            @show i
        end
        # for all potential parents Pᵢ ∈ 𝒫s

        #for
    end
    return Is
end
