export GraphAlgorithm
export infer_graph

"""
    GraphAlgorithm

The supertype of all causal graph inference algorithms.

## Concrete implementations

- [`OCE`](@ref). The optimal causation entropy algorithm for time series graphs.
- [`PC`](@ref).
"""
abstract type GraphAlgorithm end

"""
    infer_graph(algorithm::GraphAlgorithm, x) → g

Infer graph from input data `x` using the given `algorithm`.

Returns `g`, whose type depends on `algorithm`.
"""
function infer_graph end
