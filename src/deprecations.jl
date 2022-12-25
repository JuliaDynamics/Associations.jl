
function crossmap(x::AbstractVector, y::AbstractVector, args...; kwargs...)
    error("""
    The cross-mapping API has been completely re-written for CausalityTools v2.0.

    This update fixes an indexing bug for an internal cross-map function, and makes it
    easier to apply different sampling approaches and run ensemble analyses.

    As part of this rewrite, the syntax `crossmap(x::AbstractVector, y::AbstractVector, d, τ)`
    has been phased out without deprecation. The new syntax is, for target `t` (the time
    series being predicted) and source time series `s` (whose embedding is used for
    neighbor searches to predict `t`), is:
    - `crossmap(ConvergentCrossMapping(; d, τ, ...), t, s)`

    `crossmap_bootstrap` has been replaced with `crossmap(::CrossmapEnsemble, x, y)`.

    See the online documentation for more information.
    """)
end

function pai(x::AbstractVector, y::AbstractVector, args...; kwargs...)
    error("""
    The cross-mapping API has been completely re-written for CausalityTools v2.0.

    This update fixes an indexing bug for an internal cross-map function, and makes it
    easier to apply different sampling approaches and run ensemble analyses.

    As part of this rewrite, the syntax `pai(x::AbstractVector, y::AbstractVector, d, τ)`
    has been phased out without deprecation. The new syntax is, for target `t` (the time
    series being predicted) and source time series `s` (which is included as a non-lagged
    variable in the embedding of `t`), is:
    - `crossmap(PairwiseAsymmetricInference(; d, τ, ...), t, s)`

    `crossmap_bootstrap` has been replaced with `crossmap(::CrossmapEnsemble, x, y)`.

    See the online documentation for more information.
    """)
end
