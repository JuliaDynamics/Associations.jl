using DataStructures: SortedDict

include("pcmci_utils.jl")

export PCMCI
export PCMCIParent

"""
    PCMCI <: GraphAlgorithm
    PCMCI(; τmax::Int = 5, pmax::Int = -1, qmax::Int = 1, pX::Int = 3, α = 0.05,
        fdr_adjust = true,
        utest = SurrogateAssociationTest(KSG2(MIShannon(), k = 10), nshuffles = 19),
        ctest = LocalPermutationTest(MesnerShalizi(CMIShannon(), k=10, w=0), nshuffles = 19),
    )

The PCMCI causal network inference algorithm [Runge2019](@cite). 

# Keyword arguments

The PCMCI algorithm is here implemented with compatibility with *any* of the pairwise and 
conditional independence tests implemented in Associations.jl.

- **`utest`**: The pairwise/unconditional association test. Defaults to a [`SurrogateAssociationTest`](@ref)
    with a mutual information estimator, but can be any [`IndependenceTest`](@ref).
- **`ctest`**: The conditional association test. Defaults to a [`SurrogateAssociationTest`](@ref)
    with a conditional mutual information estimator, but can be any conditional [`IndependenceTest`](@ref).
- **`α`**: Significance threshold. 
- **`τmax`**: The maximum lag to consider. 
- **`pmax`**: The maximum dimension of the condition. Use `-1` (the default) for no
    restriction, corresponding to the paper's `pmax = Nτmax`, where `N` is the number of
    variables (`dimension(input_dataset)`) and `τmax` is the maximum lag.
- **`qmax`**: The maximum number of combinations.
- **`pX`**: The maximum number of the driver variable's strongest parents to condition on in the MCI stage 
    (last term of Eq. 3 in [Runge2019](@cite)). Conditioning on these accounts for autocorrelation. Too low 
    risks `pX` inflated false positives, too high lowers test power and increases dimensionality (and hence 
    increases computation time).
- **`fdr_adjust`**: Adjust p-values for all links using the False Discovery Rate (FDR) approach?
- **`use_abs_utest`**: Whether to rank parents by the absolute value `|I|` of the *pairwise*
    (`utest`) test statistic (`true`), or by the raw value (`false`). Defaults to `true`.
- **`use_abs_ctest`**: Whether to rank parents by the absolute value `|I|` of the *conditional*
    (`ctest`) test statistic (`true`), or by the raw value (`false`). Defaults to `true`.

    Using `|I|` (the paper's convention) is meaningful for signed measures like partial
    correlation, where `I = -0.8` is a *stronger* dependence than `I = +0.3`. For unsigned
    measures like MI/CMI, some estimators may yield negative values despite a theoretical
    minimum of `0`; in these cases, set the corresponding flag to `false` to compare the raw
    values directly. The flags are separate because `utest` and `ctest` may use different
    measures with different sign behaviour.

## Used with 

- [`infer_graph`](@ref). The input variables must be given as a `StateSpaceSet` or a 
    `Vector{Vector{<:Real}}`.
## Returns

When used with [`infer_graph`](@ref), it returns a [`PCMCIResult`](@ref), where
`p.parents[i]` are the parents for the variable `Xⁱₜ`. This result can be converted to a
`SimpleDiGraph` from Graphs.jl by calling `SimpleDiGraph(res)`.
"""
Base.@kwdef struct PCMCI{U,C,T} <: GraphAlgorithm
    utest::U = SurrogateAssociationTest(KSG2(MIShannon(), k=10, w=0), nshuffles=19)
    ctest::C = LocalPermutationTest(MesnerShalizi(CMIShannon(), k=10, w=0), nshuffles=19)
    τmax::T = 5
    pmax::Int = -1 # maximum condition dimension; -1 = unrestricted (N*τmax)
    qmax::Int = 1 # maximum number of combinations
    pX::Int = 3
    α = 0.05
    fdr_adjust = true
    use_abs_utest::Bool = true
    use_abs_ctest::Bool = true
end

function infer_graph(alg::PCMCI, x::AbstractStateSpaceSet; verbose=false)
    D = dimension(x)
    verbose && print_pcmci_info()

    # Find the parents of each variable in the dataset (algorithm S1)
    verbose && printstyled("˧ PC₁ condition-selection stage...\n"; color=SYMBOL_COLOR)
    𝒫 = Vector{Vector{Tuple{Int,Int}}}(undef, D)
    for i in 1:D
        𝒫[i] = first(select_parents(alg, x, i; verbose))[(i, 0)]
        verbose && print_preliminary_parents(i, 𝒫[i])
    end

    # Perform the MCI routine on the variables inferred in the first step (algorithm S2)
    verbose && printstyled("˧ MCI causal-discovery stage...\n"; color=SYMBOL_COLOR)
    if verbose && alg.fdr_adjust
        printstyled("  (displayed p is the raw MCI p-value; final links use FDR-adjusted p, shown below)\n";
            color=SYMBOL_COLOR)
    end
    res = mci_causal_discovery(alg, x, 𝒫; verbose)
    if verbose
        printstyled("˧ Inferred parents"; color=SYMBOL_COLOR)
        alg.fdr_adjust && printstyled(" (p-values FDR-adjusted)"; color=SYMBOL_COLOR)
        printstyled(":\n"; color=SYMBOL_COLOR)
        show(stdout, MIME"text/plain"(), res)
        println() # separate the verbose progress block from any caller/REPL echo of `res`
        print_pcmci_summary(alg, res, D)
    end
    return res
end

# Print a one-line summary of the inferred graph: the number of links found (including
# auto-dependencies) and the key parameters used.
function print_pcmci_summary(alg::PCMCI, res::PCMCIResult, D::Int)
    nlinks = sum(length, res.parents)
    printstyled("˧ Found $nlinks link$(nlinks == 1 ? "" : "s") (incl. auto-dependencies) among $D variables "; color=SYMBOL_COLOR)
    printstyled("(τmax=$(alg.τmax), α=$(alg.α)$(alg.fdr_adjust ? ", FDR-adjusted" : "")).\n"; color=SYMBOL_COLOR)
end

# Print a legend explaining the notation used in the verbose PCMCI output. This mirrors
# the info message printed at the start of the OCE algorithm (see `OCEInfoMessage`).
function print_pcmci_info()
    printstyled("Inferring causal graph using PCMCI\n"; bold=true)
    printstyled("Notation:\n"; underline=true, color=:default)
    printstyled("  a ⫫ b | c  := `a` is conditionally independent of `b`, given `c`\n";
        color=:default)
    printstyled("  a !⫫ b | c := `a` is conditionally dependent of `b`, given `c`\n";
        color=:default)

    # Target variable.
    print_lagged("* xᵢ", "τ"; color=TARGET_COLOR)
    printstyled(" := target variable at lag "; color=:default)
    printstyled("τ\n"; color=LAG_COLOR)

    # Driver (candidate parent) variable.
    print_lagged("* xⱼ", "τ"; color=SOURCE_COLOR)
    printstyled(" := driver (candidate parent) variable at lag "; color=:default)
    printstyled("τ\n"; color=LAG_COLOR)

    # Conditioning set.
    print_lagged("* 𝒫ᵢ", "τ"; color=CONDITIONAL_COLOR)
    printstyled(" := conditioning set (selected parents) of "; color=:default)
    print_lagged("xᵢ", "τ"; color=TARGET_COLOR)
    print("\n")
end

# Print the preliminary parents `𝒫` (a vector of `(index, lag)` tuples) selected for
# variable `Xⁱ` in the PC₁ condition-selection stage.
function print_preliminary_parents(i::Int, 𝒫::Vector{Tuple{Int,Int}})
    print_lagged("  " * add_subscript("x", i), 0; color=TARGET_COLOR)
    printstyled(" ← "; color=SYMBOL_COLOR)
    print_condset(𝒫)
    println()
end

# Print a conditioning set given as a vector of `(index, lag)` tuples, formatted as
# `{xⱼ(τ), ...}`. An empty set is printed as `∅`.
function print_condset(𝒮::Vector{Tuple{Int,Int}})
    if isempty(𝒮)
        printstyled("∅"; color=CONDITIONAL_COLOR)
        return
    end
    n = length(𝒮)
    printstyled("{"; color=CONDITIONAL_COLOR)
    for (r, (j, τ)) in enumerate(𝒮)
        print_lagged(add_subscript("x", j), τ; color=CONDITIONAL_COLOR)
        r < n && printstyled(", "; color=SYMBOL_COLOR)
    end
    printstyled("}"; color=CONDITIONAL_COLOR)
end

# Print the header for one PC₁ iteration of condition dimension `p`, over `ncandidates`
# remaining candidate parents. `p = 0` tests unconditional (pairwise) associations; `p > 0`
# conditions on `p` of the currently strongest candidates.
function print_pc1_iteration_header(p::Int, ncandidates::Int)
    printstyled("  ┌ p=$p: "; color=SYMBOL_COLOR)
    if p == 0
        printstyled("unconditional tests"; color=SYMBOL_COLOR)
    else
        printstyled("conditioning on $p variable$(p == 1 ? "" : "s")"; color=SYMBOL_COLOR)
    end
    printstyled(" ($ncandidates candidate$(ncandidates == 1 ? "" : "s"))\n"; color=SYMBOL_COLOR)
end

# Print the footer for one PC₁ iteration: the candidates `𝒫` that survived, sorted from
# strongest to weakest dependence.
function print_pc1_iteration_footer(𝒫::Vector{Tuple{Int,Int}})
    n = length(𝒫)
    printstyled("  └ $n candidate$(n == 1 ? " remains" : "s remain"): "; color=SYMBOL_COLOR)
    print_condset(𝒫)
    println()
end

# Print a single PC₁ removal decision: the candidate parent xⱼ(τ) was found conditionally
# independent of the target xᵢ(0) given the tested subset `𝒮`, and is thus removed from the
# preliminary parent set 𝒫(xᵢ).
function print_pc1_removal(i::Int, j::Int, τ::Int, 𝒮::Vector{Tuple{Int,Int}})
    print_lagged("  │   " * add_subscript("x", i), 0; color=TARGET_COLOR)
    printstyled(" ⫫ "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), τ; color=SOURCE_COLOR)
    printstyled(" | "; color=SYMBOL_COLOR)
    print_condset(𝒮)
    printstyled(" → removing "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), τ; color=SOURCE_COLOR)
    printstyled("\n"; color=SYMBOL_COLOR)
end

# Print a single MCI test decision for the potential link xⁱ(-τ) → xʲ(0), conditioned on
# the set `Z`, together with the p-value and test-statistic value.
function print_mci_test(alg::PCMCI, i::Int, j::Int, τ::Int, Z::Vector{Tuple{Int,Int}}, pval, I)
    significant = pval < alg.α
    depsymb = significant ? " !⫫ " : " ⫫ "
    print_lagged("    " * add_subscript("x", i), -τ; color=SOURCE_COLOR)
    printstyled(depsymb; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), 0; color=TARGET_COLOR)
    printstyled(" | "; color=SYMBOL_COLOR)
    print_condset(Z)
    printstyled(" [p=$(round(pval; digits=4)), I=$(round(I; digits=4))]"; color=SYMBOL_COLOR)
    significant && printstyled(" (candidate link)"; color=SOURCE_COLOR)
    printstyled("\n"; color=SYMBOL_COLOR)
end

# Print a consolidated summary of all candidate links (raw MCI p-value < α) gathered during
# the MCI stage, *before* any FDR adjustment. The final graph is a subset of these: FDR only
# raises p-values, so some candidates may be dropped when `fdr_adjust = true`.
function print_mci_candidate_summary(alg::PCMCI, p_values, test_statistics)
    candidates = [key for key in keys(p_values) if p_values[key] < alg.α]
    # Sort by target variable first (key[2]), then source (key[1]), then lag (key[3]).
    sort!(candidates, by=key -> (key[2], key[1], key[3]))
    n = length(candidates)
    printstyled("˧ $n candidate link$(n == 1 ? "" : "s") (raw p < α"; color=SYMBOL_COLOR)
    alg.fdr_adjust && printstyled(", before FDR adjustment"; color=SYMBOL_COLOR)
    printstyled("):\n"; color=SYMBOL_COLOR)
    for key in candidates
        i, j, τ = key  # driver `i` at lag `τ` (≤ 0), target `j` at lag 0
        print_lagged("    " * add_subscript("x", i), τ; color=SOURCE_COLOR)
        printstyled(" → "; color=SYMBOL_COLOR)
        print_lagged(add_subscript("x", j), 0; color=TARGET_COLOR)
        printstyled(" [p=$(round(p_values[key]; digits=4)), I=$(round(test_statistics[key]; digits=4))]\n";
            color=SYMBOL_COLOR)
    end
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
    # Initialize dictionary to store p-values and test statistic values for causal links
    p_values = SortedDict{Tuple{Int,Int,Int},Float64}()
    test_statistics = SortedDict{Tuple{Int,Int,Int},Float64}()

    # Iterate over all potential links between X^i and X^j, considering different time lags
    D = dimension(X)


    for i = 1:D # driver variable Xⁱ (tested at lag t-τ below)
        𝒫xₜⁱ = 𝒫[i]
        for j = 1:D # target variable Xʲ (tested at lag t below)
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
                    𝒫pxXⁱt₋τ = [(i, τᵢ - τ) for (i, τᵢ) in 𝒫pxXⁱt]

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
                    verbose && print_mci_test(alg, i, j, τ, Z, pval, I)
                end
            end
        end
    end

    verbose && print_mci_candidate_summary(alg, p_values, test_statistics)

    if alg.fdr_adjust
        p_values_adjusted = fdr_adjust(alg, p_values)
    else
        p_values_adjusted = p_values
    end

    # Format results.
    D = dimension(X)
    parents = Vector{Vector{PCMCIParent}}(undef, D)
    for i in 1:D
        parents[i] = [
            PCMCIParent(key[1], key[3], p_values_adjusted[key], test_statistics[key]) for key in keys(p_values_adjusted) if key[2] == i && p_values_adjusted[key] < alg.α
        ]
    end
    return PCMCIResult(parents)
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
                if pval > alg.α
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

"""
    lexicographical_subsets(set::Vector{Tuple{Int, Int}}, cardinality::Int; 
        exclude::Tuple{Int, Int} = nothing)

Generating lexicographical subsets with `cardinality` from the given `set` of 
integer tuples of the form `(i, τ)`, where `i` is the index of the variable
(in the original `StateSpaceSet` given to `infer_graph`) , and `τ` is the 
lag. If a tuple is given as `exclude`, then this `(index, lag)` is excluded
from the vector of returned subsets.
"""
function lexicographical_subsets(set::Vector{Tuple{Int,Int}}, cardinality::Int; exclude::Tuple{Int,Int}=nothing)
    filtered_set = isnothing(exclude) ? set : filter(x -> x != exclude, set)
    return [subset for subset in combinations(filtered_set, cardinality)]
end


"""
    fdr_adjust(alg::PCMCI, p_values::SortedDict{Tuple{Int, Int, Int}, Float64})

Adjust `p_values` to control false discovery rate using the Benjamini-Hochberg (1995)
procedure. 

## Arguments
- **`p_values`*:  Dictionary of p-values with keys representing (i, j, τ)
- **`α`**: Significance level.
"""
function fdr_adjust(alg::PCMCI, p_values::SortedDict{Tuple{Int,Int,Int},Float64})
    # Extract p-values and sort them in ascending order
    sorted_keys = collect(keys(p_values))
    sorted_p_values = collect(values(p_values))
    n_tests = length(sorted_p_values)

    # Sort the p-values and keep track of the original indices
    sorted_indices = sortperm(sorted_p_values)
    sorted_p_values = sorted_p_values[sorted_indices]
    sorted_keys = sorted_keys[sorted_indices]

    # Apply the BH procedure to compute adjusted p-values
    adjusted_p_values = Vector{Float64}(undef, n_tests)
    for i in 1:n_tests
        rank = i
        # Compute the BH critical value
        adjusted_p_values[i] = min(1.0, sorted_p_values[i] * n_tests / rank)
    end

    # Ensure adjusted p-values are non-decreasing
    for i in (n_tests-1):-1:1
        adjusted_p_values[i] = min(adjusted_p_values[i], adjusted_p_values[i+1])
    end

    # Return a dictionary with similar structure, but with the adjusted p-values
    adjusted_p_values_dict = SortedDict{Tuple{Int,Int,Int},Float64}()
    for i in 1:n_tests
        adjusted_p_values_dict[sorted_keys[i]] = adjusted_p_values[i]
    end

    return adjusted_p_values_dict
end
