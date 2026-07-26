using DataStructures: SortedDict

include("pcmci_types.jl")
include("pcmci_utils.jl")
include("pcmci_printing.jl")

function infer_graph(alg::PCMCI, x::AbstractStateSpaceSet; verbose=false)
    D = dimension(x)
    if verbose
        print_pcmci_info()
        print_pcmci_parameters(alg)
    end

    # Algorithm S1: find the parents of each variable in the dataset.
    verbose && printstyled("˧ PC₁ condition-selection stage...\n"; color=SYMBOL_COLOR)
    𝒫 = Vector{Vector{Tuple{Int,Int}}}(undef, D)
    for i in 1:D
        𝒫[i] = first(select_parents(alg, x, i; verbose))[(i, 0)]
        verbose && print_preliminary_parents(i, 𝒫[i])
    end

    # Algorithm S2: Perform the MCI routine on the variables inferred in the first step
    res = mci_causal_discovery(alg, x, 𝒫; verbose)
    return res
end
function infer_graph(alg::PCMCI, x::Vector{Vector{T}}; kwargs...) where T
    return infer_graph(alg, StateSpaceSet(x...); kwargs...)
end

"""
    mci_causal_discovery(alg::PCMCI, X::AbstractStateSpaceSet, 𝒫)

Algorithm S2 in the original paper.

## Arguments

- **`alg`**: The algorithm with parameters for the analysis.
- **`X`**: `StateSpaceSet` containing the multivariate time series data.
- **`𝒫`**: Preliminary parents obtained from Algorithm S1 for all variables in `X`.
"""
function mci_causal_discovery(alg::PCMCI, X::AbstractStateSpaceSet, 𝒫; verbose=false)
    verbose && print_mci_stage_header(alg)

    # Initialize dictionary to store p-values and test statistic values for causal links
    p_values = SortedDict{Tuple{Int,Int,Int},Float64}()
    test_statistics = SortedDict{Tuple{Int,Int,Int},Float64}()

    # Iterate over all potential links between X^i and X^j, considering different time lags
    D = dimension(X)

    for j = 1:D # target variable Xʲ (tested at lag t below)
        survivors = verbose ? Tuple{Int,Int}[] : nothing
        verbose && print_mci_target_header(alg, j, D)
        for i = 1:D # driver variable Xⁱ (tested at lag t-τ below)
            𝒫xₜⁱ = 𝒫[i]
            for τ in 0:alg.τmax
                # A link Xⁱt₋τ → Xʲₜ exists iff Xⁱt₋τ !⫫ Xʲₜ | 𝒫(Xtʲ) \ {Xⁱt₋τ}, 𝒫px(Xⁱt₋τ)
                if τ > 0 || i != j  # Avoid self-links for contemporaneous relationships
                    # ----------------------------------------------------------------
                    # Define conditioning set
                    # ----------------------------------------------------------------
                    𝒫xₜʲ = filter(x -> x != (i, -τ), 𝒫[j])

                    # Get first `n` relevant parent variables of Xⁱt₋τ, limited by `pX::Int`
                    n = min(alg.pX, length(𝒫xₜⁱ))
                    𝒫pxXⁱt = first(𝒫xₜⁱ, n)

                    # Time-shift these variables by the outer (loop-based) τ
                    𝒫pxXⁱt₋τ = [(k, τᵢ - τ) for (k, τᵢ) in 𝒫pxXⁱt]

                    # Finally combine the variables into the conditioning set Z
                    Z::Vector{Tuple{Int,Int}} = union(𝒫xₜʲ, 𝒫pxXⁱt₋τ)

                    # Add conditioning variables from the parents set to `ts` and `js`
                    ts = [0, -τ]  # Embedding time delays: current time 0 and lag -τ
                    js = [j, i]  # Indices of the variables involved: target j and lagged driver i
                    cond_τs = [τ for (i, τ) in Z]
                    cond_idxs = [i for (i, τ) in Z]
                    append!(ts, cond_τs)
                    append!(js, cond_idxs)

                    # Generate the embedding and test (conditional) independence between
                    # X^i_{t-τ} and X^j_t.
                    embedding = genembed(X, ts, js)
                    pval, I = independence_test_on_subset(alg, embedding)
                    p_values[(i, j, -τ)] = pval
                    test_statistics[(i, j, -τ)] = I
                    if verbose
                        print_mci_test(alg, i, j, τ, Z, 𝒫xₜʲ, pval, I; prefix="  │   ")
                        pval < alg.α && push!(survivors, (i, -τ))
                    end
                end
            end
        end
        verbose && print_mci_target_footer(j, survivors)
    end

    verbose && print_mci_candidate_summary(alg, p_values, test_statistics)

    if alg.fdr_adjust
        p_values_adjusted = fdr_adjust(alg, p_values)
    else
        p_values_adjusted = p_values
    end

    # Format results: one `PCMCISelectedParents` per variable (`res[i]` are the parents of
    # `xᵢ(0)`), matching the per-variable return convention of the OCE algorithm.
    D = dimension(X)
    res = Vector{PCMCISelectedParents}(undef, D)
    for i in 1:D
        links = [
            PCMCIParent(key[1], key[3], p_values_adjusted[key], test_statistics[key]) for key in keys(p_values_adjusted) if key[2] == i && p_values_adjusted[key] < alg.α
        ]
        res[i] = PCMCISelectedParents(i, links)
    end
    verbose && print_inferred_parents(alg, res, D)
    return res
end

"""
    select_parents(alg::PCMCI, x, i::Int)

Select the parents of the `i`-th variable in `x` according to the `PCMCI` algorithm.
This is algorithm S1 in their supplementary material.
"""
function select_parents(alg::PCMCI, x::AbstractStateSpaceSet, i::Int; verbose=false)
    if verbose
        printstyled("  Selecting parents for "; color=SYMBOL_COLOR)
        print_lagged(add_subscript("x", i), 0; color=TARGET_COLOR)
        printstyled("...\n"; color=SYMBOL_COLOR)
    end

    # Preliminary parents for the `j`-th variable at lag `t-0`, limiting the
    # maximum lag to `alg.τmax`. This is a vector of integer tuple of the form `(i, τ)`.
    # Many of the elements of `𝒫ₜʲ` will be gradually eliminated in the loop below.
    𝒫ₜʲ = initialize_parents_for_variable(i, dimension(x), alg.τmax)

    # Dictionary of test statistics Imin(Xⁱ_{t-tau} - Xʲ_{t}) = ∞ for all parents in 𝒫ₜʲ
    Imin = initialize_test_statistics_for_variable(i, dimension(x), alg.τmax)

    # `pmax = -1` means unrestricted, i.e. the paper's `pmax = Nτmax`.
    pmax = alg.pmax < 0 ? dimension(x) * alg.τmax : alg.pmax

    for p = 0:pmax
        if length(𝒫ₜʲ[(i, 0)]) - 1 < p
            break  # This will exit the entire loop
        end
        verbose && print_pc1_iteration_header(p, length(𝒫ₜʲ[(i, 0)]))
        # Iterate over all potential parents.
        marked_for_removal = Vector{Tuple{Int,Int}}(undef, 0)
        for (j, τ) in 𝒫ₜʲ[(i, 0)]
            q = -1
            # All lexicographically ordered subsets of cardinality `p` from 𝒫ₜʲ
            for variable_subset in lexicographical_subsets(𝒫ₜʲ[(i, 0)], p, exclude=(j, τ))
                q = q + 1
                if q >= alg.qmax
                    break
                end

                if p == 0
                    # Embedding: 
                    # - the first column is Xʲₜ
                    # - the second column is Xʲτ
                    embedding = genembed(x, (0, τ), (i, j))
                else
                    varinds = [j for (j, τ) in variable_subset]
                    varτs = [τ for (j, τ) in variable_subset]

                    # Embedding: 
                    # - the first column is Xʲₜ
                    # - the second column is Xʲτ
                    # - the third-to-end columns are the variables in variable_subset
                    embedding = genembed(x, (0, τ, varτs...), (i, j, varinds...))
                end
                pval, I = independence_test_on_subset(alg, embedding)

                # Track minimum dependence strength across the tested conditioning
                # subsets (Algorithm S1, lines 16-17).
                #
                # With `use_abs = true` (the paper's convention) strength is measured by |I|.
                # This is the correct notion for a signed measure (e.g. partial correlation),
                # where I = -0.8 is a *stronger* dependence than I = +0.3. For e.g. MI/CMI
                # estimators that may yield negative numbers despite the theoretical minimum
                # actually being 0, this is not the correct behaviour: then smaller negative
                # numbers are actually smaller.
                #
                # `p == 0` uses the pairwise `utest`; `p > 0` uses the conditional `ctest`,
                # so we pick the corresponding `use_abs` flag.
                use_abs = p == 0 ? alg.use_abs_utest : alg.use_abs_ctest
                Iₘ = use_abs ? abs(I) : I
                if Iₘ < Imin[(i, 0)][(j, τ)]
                    Imin[(i, 0)][(j, τ)] = Iₘ
                end
                if pval > alg.αPC
                    push!(marked_for_removal, (j, τ))
                    verbose && print_pc1_removal(i, j, τ, variable_subset)
                    break
                end
            end
        end
        remove_marked_potential_parents!(𝒫ₜʲ, marked_for_removal, i)

        # Algorithm S1, line 22 / Materials and Methods: after every iteration, sort the
        # remaining preliminary parents by their `Imin` (minimum |I| across tested subsets),
        # strongest dependence first.
        sort!(𝒫ₜʲ[(i, 0)], by=parent -> Imin[(i, 0)][parent], rev=true)
        verbose && print_pc1_iteration_footer(𝒫ₜʲ[(i, 0)])
    end
    return 𝒫ₜʲ, Imin
end

# We need to pass a different number of arguments to `independence` depending on whether 
# the embedding is two-dimension or higher-dimensional.
function independence_test_on_subset(alg::PCMCI, embedding::AbstractStateSpaceSet{D,T}) where {D,T}
    if D == 2
        Xʲₜ = embedding[:, 1] |> StateSpaceSet
        Xʲτ = embedding[:, 2] |> StateSpaceSet
        test_result = independence(alg.utest, Xʲₜ, Xʲτ)
    else
        Xʲₜ = embedding[:, 1] |> StateSpaceSet
        Xʲτ = embedding[:, 2] |> StateSpaceSet
        if D == 3
            X_variable_subset = embedding[:, 3:D] |> StateSpaceSet
        else
            X_variable_subset = embedding[:, 3:D]
        end
        test_result = independence(alg.ctest, Xʲₜ, Xʲτ, X_variable_subset)
    end

    pval = pvalue(test_result)
    I = test_statistic(test_result)

    return pval, I
end

function remove_marked_potential_parents!(𝒫ₜʲ, marked_for_removal::Vector{Tuple{Int,Int}}, j)
    𝒫ₜʲ[(j, 0)] = filter(x -> x ∉ marked_for_removal, 𝒫ₜʲ[(j, 0)])
    return 𝒫ₜʲ
end

function initialize_parents_for_variable(j::Int, num_variables::Int, max_lag::Int)
    parents = SortedDict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    # Initialize only the parents for X_t^j (not lagged versions of X^j itself)
    parents[(j, 0)] = [
        (i, -tau) for i in 1:num_variables, tau in 1:max_lag
                      if i != j || tau != 0
    ]
    return parents
end

function initialize_test_statistics_for_variable(j::Int, num_variables::Int, max_lag::Int)
    test_statistics = SortedDict{Tuple{Int,Int},Dict{Tuple{Int,Int},Float64}}()
    test_statistics[(j, 0)] = Dict{Tuple{Int,Int},Float64}()
    # Initialize each test statistic value to Inf for all potential parents of X_t^j
    for parent in [(i, -tau) for i in 1:num_variables, tau in 1:max_lag if i != j || tau != 0]
        test_statistics[(j, 0)][parent] = Inf
    end
    return test_statistics
end

