# Cross entropy

Cross entropies appear in many forms in the literature, depending on what type of entropy
they are based on. In CausalityTools.jl, we indicate the different types of cross entropies
by passing an [`CrossEntropy`](@ref) instance as the first argument to 
[`entropy_cross`](@ref).

For example, to compute the Shannon cross entropy, do
`entropy_cross(CrossEntropyShannon(), est, x)`, and to compute Rényi cross entropy, do
`entropy_cross(CrossEntropyRenyi(; q), est, x)`.

## Implementations

Any of the following estimators can be used with [`entropy_cross`](@ref).

| Estimator                  | Principle         | Input data                 | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| -------------------------- | ----------------- | -------------------------- | :---------------: | :-------------: | :---------------: |
| [`BulinskiDimitrov`](@ref) | Nearest neighbors | Equal-dimension `Dataset`s |     Discrete      |        x        |         x         |

## API

```@docs
entropy_cross
CrossEntropy
CrossEntropyDefinition
CrossDifferentialEntropyEstimator
```

### Shannon cross-entropy

```@docs
CrossEntropyShannon
BulinskiDimitrov
```

### Rényi cross-entropy

```@docs
CrossEntropyRenyi
Thierrin
```
