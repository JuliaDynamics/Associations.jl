# In an MCI conditioning set (Eq. 3) the conditions have two origins. The target's parents
# P̂(Xʲ)\{driver} use CONDITIONAL_COLOR (as in PC₁), while the driver's parents P̂_pX(Xⁱ),
# added for autocorrelation control, use CONDITIONAL_DRIVER_COLOR — tinted like the driver
# (SOURCE_COLOR) so they read as "conditions contributed by the driver".
const CONDITIONAL_DRIVER_COLOR = SOURCE_COLOR

# Print a one-line summary of the inferred graph: the number of links found (including
# auto-dependencies) and the key parameters used.
function print_pcmci_summary(alg::PCMCI, res::AbstractVector{<:PCMCISelectedParents}, D::Int)
    nlinks = sum(nparents, res)
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

    # Conditioning-set colour key for the MCI stage (Eq. 3).
    printstyled("* MCI conditioning set combines "; color=:default)
    printstyled("target's parents"; color=CONDITIONAL_COLOR)
    printstyled(" and the "; color=:default)
    printstyled("driver's parents"; color=CONDITIONAL_DRIVER_COLOR)
    printstyled(" (added for autocorrelation control)\n"; color=:default)
end

# Compact label for an independence test: outer wrapper name plus its first type parameter,
# which is always the wrapped measure/estimator, e.g. `SurrogateAssociationTest{ChatterjeeCorrelation}`.
function test_label(test)
    T = typeof(test)
    outer = nameof(T)
    params = T.parameters
    isempty(params) && return string(outer)
    inner = first(params)
    inner_name = inner isa Type ? nameof(inner) : nameof(typeof(inner))
    return string(outer, "{", inner_name, "}")
end

# Print the algorithm parameters used for this run, so the verbose output is self-contained.
function print_pcmci_parameters(alg::PCMCI)
    printstyled("Parameters:\n"; underline=true, color=:default)
    # Print one parameter line: grey name/value, grey parenthetical description.
    function param(name, value, desc)
        printstyled("  $(rpad(name, 13)) = "; color=SYMBOL_COLOR)
        printstyled(value; color=:default)
        printstyled("  ($desc)\n"; color=SYMBOL_COLOR)
    end
    pmax_str = alg.pmax < 0 ? "-1 (unrestricted, i.e. N·τmax)" : string(alg.pmax)
    param("utest", test_label(alg.utest), "unconditional test; PC₁ p=0 & MCI unconditional links")
    param("ctest", test_label(alg.ctest), "conditional test; PC₁ p>0 & MCI conditional links")
    param("τmax", alg.τmax, "maximum time lag considered")
    param("pmax", pmax_str, "max PC₁ conditioning-set dimension")
    param("qmax", alg.qmax, "max conditioning subsets tested per candidate in PC₁")
    param("pX", alg.pX, "driver's strongest parents added to MCI conditioning set")
    param("α", alg.α, "significance level")
    param("fdr_adjust", alg.fdr_adjust, "Benjamini-Hochberg FDR adjustment of MCI p-values")
    param("use_abs_utest", alg.use_abs_utest, "rank PC₁ p=0 dependence by |I|")
    param("use_abs_ctest", alg.use_abs_ctest, "rank PC₁ p>0 dependence by |I|")
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
# preliminary parent set 𝒫(xᵢ). The candidate is printed first (`xⱼ(τ) ⫫ xᵢ(0)`), matching
# the driver-first order of the paper's Eqs. 3-4 and the MCI stage.
function print_pc1_removal(i::Int, j::Int, τ::Int, 𝒮::Vector{Tuple{Int,Int}})
    print_lagged("  │   " * add_subscript("x", j), τ; color=SOURCE_COLOR)
    printstyled(" ⫫ "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", i), 0; color=TARGET_COLOR)
    printstyled(" | "; color=SYMBOL_COLOR)
    print_condset(𝒮)
    printstyled(" → removing "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), τ; color=SOURCE_COLOR)
    printstyled(" as candidate parent of "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", i), 0; color=TARGET_COLOR)
    printstyled("\n"; color=SYMBOL_COLOR)
end

# Print a single MCI test decision for the potential link xⁱ(-τ) → xʲ(0), conditioned on
# the set `Z`, together with the p-value and test-statistic value.
# Print an MCI conditioning set `Z`, colouring each member by its origin (Eq. 3): members
# that are among the target's parents `Zt` use CONDITIONAL_COLOR (as in PC₁), while the
# remaining members — the driver's pX parents added for autocorrelation control — use
# CONDITIONAL_DRIVER_COLOR. A member present in both origins is shown as a target parent.
function print_mci_condset(Z::Vector{Tuple{Int,Int}}, Zt::Vector{Tuple{Int,Int}})
    if isempty(Z)
        printstyled("∅"; color=CONDITIONAL_COLOR)
        return
    end
    n = length(Z)
    printstyled("{"; color=CONDITIONAL_COLOR)
    for (r, (j, τ)) in enumerate(Z)
        c = (j, τ) in Zt ? CONDITIONAL_COLOR : CONDITIONAL_DRIVER_COLOR
        print_lagged(add_subscript("x", j), τ; color=c)
        r < n && printstyled(", "; color=SYMBOL_COLOR)
    end
    printstyled("}"; color=CONDITIONAL_COLOR)
end

function print_mci_test(alg::PCMCI, i::Int, j::Int, τ::Int, Z::Vector{Tuple{Int,Int}}, Zt::Vector{Tuple{Int,Int}}, pval, I; prefix::AbstractString="    ")
    significant = pval < alg.α
    depsymb = significant ? " !⫫ " : " ⫫ "
    print_lagged(prefix * add_subscript("x", i), -τ; color=SOURCE_COLOR)
    printstyled(depsymb; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), 0; color=TARGET_COLOR)
    printstyled(" | "; color=SYMBOL_COLOR)
    print_mci_condset(Z, Zt)
    printstyled(" [p=$(round(pval; digits=4)), I=$(round(I; digits=4))]"; color=SYMBOL_COLOR)
    # State the decision explicitly: is the driver xᵢ(-τ) kept (dependence survives
    # conditioning → candidate link) or discarded (conditionally independent) as a
    # parent of the target xⱼ(0)?
    printstyled(significant ? " → keeping " : " → discarding "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", i), -τ; color=SOURCE_COLOR)
    printstyled(significant ? " as parent of " : " as candidate parent of "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), 0; color=TARGET_COLOR)
    printstyled("\n"; color=SYMBOL_COLOR)
end

# Print the header announcing the start of the MCI causal-discovery stage (Algorithm S2). When
# FDR adjustment is enabled, note that the per-test p-values displayed below are the raw MCI
# p-values, while the final links use the FDR-adjusted p-values shown in the summary blocks.
function print_mci_stage_header(alg::PCMCI)
    printstyled("˧ MCI causal-discovery stage...\n"; color=SYMBOL_COLOR)
    if alg.fdr_adjust
        printstyled("  (displayed p is the raw MCI p-value; final links use FDR-adjusted p, shown below)\n";
            color=SYMBOL_COLOR)
    end
end

# Print the final inferred-parents block for the MCI stage: one line per variable (via the
# `AbstractSelectedParents` `show` method) followed by a one-line summary. When FDR adjustment
# is enabled, the reported p-values are the FDR-adjusted ones.
function print_inferred_parents(alg::PCMCI, res::AbstractVector{<:PCMCISelectedParents}, D::Int)
    printstyled("˧ Inferred parents"; color=SYMBOL_COLOR)
    alg.fdr_adjust && printstyled(" (p-values FDR-adjusted)"; color=SYMBOL_COLOR)
    printstyled(":\n"; color=SYMBOL_COLOR)
    for p in res
        show(stdout, MIME"text/plain"(), p)
        println()
    end
    println() # separate the verbose progress block from any caller/REPL echo of `res`
    print_pcmci_summary(alg, res, D)
end

# Open a boxed per-target group for the MCI stage, mirroring the per-variable layout of the
# PC₁ stage. The target-major MCI loop makes all tests for `xⱼ(0)` contiguous, so decisions
# are printed inline (between this header and `print_mci_target_footer`) with no buffering.
# The number of tests is fixed: every (driver, lag) pair except the self-contemporaneous
# `(j, 0)`, i.e. `D·(τmax+1) - 1`.
function print_mci_target_header(alg::PCMCI, j::Int, D::Int)
    printstyled("  Testing candidate parents of "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), 0; color=TARGET_COLOR)
    printstyled("...\n"; color=SYMBOL_COLOR)
    n = D * (alg.τmax + 1) - 1
    printstyled("  ┌ $n test$(n == 1 ? "" : "s")\n"; color=SYMBOL_COLOR)
end

# Close a per-target group, listing the drivers whose dependence survived (raw p < α). The
# final graph uses FDR-adjusted p-values (see the candidate summary and inferred-parents
# block below).
function print_mci_target_footer(j::Int, survivors::Vector{Tuple{Int,Int}})
    printstyled("  └ "; color=SYMBOL_COLOR)
    print_lagged(add_subscript("x", j), 0; color=TARGET_COLOR)
    printstyled(" ← "; color=SYMBOL_COLOR)
    print_condset(survivors)
    printstyled("   (raw p<α; final set after FDR below)\n"; color=SYMBOL_COLOR)
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

