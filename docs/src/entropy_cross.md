# Cross entropy

## API

```@docs
entropy_cross
CrossEntropyEstimator
```

## Estimators

```@docs
BulinskiDimitrov
```

## Implementations

| Estimator                  | Principle         | Input data                       | [`Shannon`](@ref) | [`Renyi`](@ref) | [`Tsallis`](@ref) |
| -------------------------- | ----------------- | -------------------------------- | :---------------: | :-------------: | :---------------: |
| [`BulinskiDimitrov`](@ref) | Nearest neighbors | `Dataset`s with equal dimensions |     Discrete      |        x        |         x         |
