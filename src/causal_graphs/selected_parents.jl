using Graphs: add_edge!
using Graphs.SimpleGraphs: SimpleDiGraph

export AbstractSelectedParents
export target, parents_js, parents_τs, nparents

"""
    AbstractSelectedParents

Supertype for the per-variable parent containers returned by graph-inference algorithms.

When an algorithm is used with [`infer_graph`](@ref), a `Vector{<:AbstractSelectedParents}` 
is returned --- one element per input variable, where `p[i]` holds the parents of the 
`i`-th variable `xᵢ(0)`.

## Interface

Subtypes provide the following fields (or overload the accessors of the same name):

- `i::Int`: index of the target variable (accessor [`target`](@ref)).
- `parents_js`: variable index of each selected parent (accessor [`parents_js`](@ref)).
- `parents_τs`: lag of each selected parent (accessor [`parents_τs`](@ref)).

Implementing this interface provides `SimpleDiGraph` conversion and console printing for
free. Concrete subtypes may carry additional per-algorithm data (e.g. per-link p-values);
override [`parent_annotation`](@ref) to include such data when printing.

## Concrete subtypes

- [`PCMCISelectedParents`](@ref) (from [`PCMCI`](@ref)).
"""
abstract type AbstractSelectedParents end

"""
    target(p::AbstractSelectedParents) → Int

The index `i` of the target variable `xᵢ(0)` whose parents are stored in `p`.
"""
target(p::AbstractSelectedParents) = p.i

"""
    parents_js(p::AbstractSelectedParents)

The variable indices of the selected parents --- one per selected parent.
"""
parents_js(p::AbstractSelectedParents) = p.parents_js

"""
    parents_τs(p::AbstractSelectedParents)

The lags of the selected parents --- one per selected parent.
"""
parents_τs(p::AbstractSelectedParents) = p.parents_τs

"""
    nparents(p::AbstractSelectedParents) → Int

The number of selected parents.
"""
nparents(p::AbstractSelectedParents) = length(parents_js(p))

# Text appended after each printed parent variable (e.g. " [p=…, I=…]"). Defaults to
# nothing; subtypes carrying per-link metadata may override.
parent_annotation(::AbstractSelectedParents, k::Int) = ""

function SimpleDiGraph(v::AbstractVector{<:AbstractSelectedParents})
    N = length(v)
    g = SimpleDiGraph(N)
    for k in 1:N
        # `v[k]` holds the parents of variable `xₖ`. Each parent is the driver `j` at lag
        # `τ`, so a link points `j → k`.
        for (j, τ) in zip(parents_js(v[k]), parents_τs(v[k]))
            if τ == 0 && j != k
                # Contemporaneous (τ = 0) links are left undirected, represented here as a
                # pair of opposing directed edges.
                add_edge!(g, j, k)
                add_edge!(g, k, j)
            elseif j != k # lagged link; skip self-loops (auto-dependencies xₖ(-τ) → xₖ(0))
                add_edge!(g, j, k)
            end
        end
    end
    return g
end

# Print the parent set of `p` as `{xⱼ(τ), ...}`, using each subtype's `parent_annotation`.
function print_condvars(p::AbstractSelectedParents)
    js, τs = parents_js(p), parents_τs(p)
    n = length(js)
    printstyled("{"; color=CONDITIONAL_COLOR)
    for r in 1:n
        print_lagged(add_subscript("x", js[r]), τs[r]; color=CONDITIONAL_COLOR)
        printstyled(parent_annotation(p, r); color=SYMBOL_COLOR)
        r < n && printstyled(", "; color=SYMBOL_COLOR)
    end
    printstyled("}"; color=CONDITIONAL_COLOR)
end

function Base.show(io::IO, ::MIME"text/plain", p::AbstractSelectedParents)
    print_lagged(add_subscript("x", target(p)), 0; color=TARGET_COLOR)
    printstyled(" ← "; color=SYMBOL_COLOR)
    print_condvars(p)
end