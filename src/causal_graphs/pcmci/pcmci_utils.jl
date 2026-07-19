export PCMCIResult

"""
    PCMCIParent(i, τ, pval, test_statistic)

The index `i` and lag `τ` of a parent node, along with the p-value and 
test statistic value for the link.
"""
struct PCMCIParent
    i::Int
    τ::Int
    pval::Real
    test_statistic::Real
end
"""
    PCMCIResult

Stores the result of a [`PCMCI`](@ref) analysis. Each parent is represented as 
a [`PCMCIParent`](@ref).

This is just a collection of `Vector{PCMCIParent}` - one per input variable.
If `r` is a `PCMCIResult`, then the parents of ``X^i_t`` are `r[i]`.

When printed in the console the `[pvalue, test_statistic]` are displayed for each
parent variable.
"""
struct PCMCIResult
    parents::Vector{Vector{PCMCIParent}}
end

function Base.show(io::IO, ::MIME"text/plain", r::PCMCIResult)
    for (i, parents::Vector{PCMCIParent}) in enumerate(r.parents)
        print_lagged(add_subscript("x", i), 0; color=TARGET_COLOR)
        printstyled(" ← "; color=SYMBOL_COLOR)
        print_condvars(parents)
        println()
    end
end

function print_condvars(parents::Vector{PCMCIParent})
    τs = [p.τ for p in parents]
    js = [p.i for p in parents]
    ps = [p.pval for p in parents]
    Is = [p.test_statistic for p in parents]

    n_selected = length(parents)
    printstyled("{", color=CONDITIONAL_COLOR)
    for r in 1:n_selected
        print_lagged(add_subscript("x", js[r]), τs[r]; color=CONDITIONAL_COLOR)
        printstyled(" [p=$(round(ps[r]; digits=4)), I=$(round(Is[r]; digits=4))]", color=SYMBOL_COLOR)
        r < n_selected && printstyled(", "; color=SYMBOL_COLOR)
    end
    printstyled("}", color=CONDITIONAL_COLOR)
end

function SimpleDiGraph(r::PCMCIResult)
    N = length(r.parents)
    g = SimpleDiGraph(N)
    for k = 1:N
        # `r.parents[k]` are the parents of variable `Xᵏ`. Each parent `p` is the driver
        # variable `p.i` at lag `p.τ`, so a lagged link points `p.i → k`.
        for p in r.parents[k]
            if p.τ == 0
                # Contemporaneous (τ = 0) links are left undirected, represented here as a
                # bidirectional edge (Runge et al. 2019, Supplementary Sect. S1). 
                add_edge!(g, p.i, k)
                add_edge!(g, k, p.i)
            elseif p.i != k # lagged link; skip self-loops (auto-dependencies xₖ(-τ) → xₖ(0))
                add_edge!(g, p.i, k)
            end
        end
    end
    return g
end
