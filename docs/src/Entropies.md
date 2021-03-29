# Generalized entropy / probabilities

The [Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl) package provides a unified interface for estimation for entropies and probabilities, using for example counting-based, nearest neighbor based, kernel density based or wavelet-based approaches.

Information theoretic causality measures are calculated using entropy estimation, so many
estimators are shared. Relevant functions are therefore reexported here for convencience.

The main **API** of Entropies.jl is contained in two functions:

* [`probabilities`](@ref) which computes probability distributions of given datasets
* [`genentropy`](@ref) which uses the output of [`probabilities`](@ref), or a set of
    pre-computed [`Probabilities`](@ref), to calculate entropies.

These functions dispatch on subtypes of [`ProbabilitiesEstimator`](@ref), which are:

```@example
using Entropies, InteractiveUtils
subtypes(ProbabilitiesEstimator)
```

## Generalized entropy

```@docs
Entropies.genentropy
```

## Probabilities

```@docs
Probabilities
probabilities
probabilities!
ProbabilitiesEstimator
```

## Fast histograms

```@docs
Entropies.binhist
```
