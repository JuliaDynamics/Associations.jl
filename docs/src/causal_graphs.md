# [Inferring causal graphs](@id causal_graphs)

Directed causal graphical models, estimated on observed data, is an incredibly
useful framework for causal inference. There exists a plethora of methods for
estimating such models.

Useful reading:

- **Pearl, J. Glymour, M., & Jewell, N. P. (2016)**. Causal inference in statistics:
    A primer. John Wiley & Sons**. An excellent introductory book, suitable for anyone
    interested, from a beginners to experts.
- **Glymour, C., Zhang, K., & Spirtes, P. (2019). Review of causal discovery methods
    based on graphical models. Frontiers in genetics, 10, 524**. The authoritative
    overview of causal discovery from graphical models. Many more methods have emerged
    since this paper, and our goal is to provide a library of these methods.

## Causal graph API

The API for inferring causal graphs is defined by:

- [`infer_graph`](@ref)
- [`GraphAlgorithm`](@ref), and its subtypes

```@docs
infer_graph
GraphAlgorithm
```

## [Optimal causation entropy](@ref)

```@docs
OCE
```

## [PC](@ref)

```@docs
PC
```
