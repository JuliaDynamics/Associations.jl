export GraphAlgorithm
export infer_graph

"""
    GraphAlgorithm

The supertype of all algorithms used to infer causal relationships from graphs.

Concrete subtypes are:

- [`PCRobust`](@ref)
"""
abstract type GraphAlgorithm end

"""
    infer_graph(algorithm::GraphAlgorithm, x) → g

Infer graph from input data `x` using the given `algorithm`.

Returns a graph `g`, whose type depends on `algorithm`.
"""
function infer_graph end
