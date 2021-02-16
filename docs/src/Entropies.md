# Entropies

The [Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl) package provides a unified interface for entropy estimation, using for example counting-based, nearest neighbor based, kernel density based or wavelet-based approaches. Some of these entropy estimators are also used to calculate mutual information and transfer entropy. Estimators relevant for 
mutual information and transfer entropy estimation are also documented here.

## Data format

Most of the code in this package assumes that your data is represented by the `Dataset`-type from [`DelayEmbeddings.jl`](https://github.com/JuliaDynamics/DelayEmbeddings.jl), where each observation is a D-dimensional data point represented by a static vector. See the [`DynamicalSystems.jl` documentation](https://juliadynamics.github.io/DynamicalSystems.jl/dev/) for more info. Univariate timeseries given as
`AbstractVector{<:Real}` also work with some estimators, but are treated differently
based on which method for probability/entropy estimation is applied.

```@docs
Dataset
```

## API

The main **API** of this package is contained in two functions:

* [`probabilities`](@ref) which computes probability distributions of given datasets
* [`genentropy`](@ref) which uses the output of [`probabilities`](@ref), or a set of
    pre-computed [`Probabilities`](@ref), to calculate entropies.

These functions dispatch on subtypes of [`ProbabilitiesEstimator`](@ref), which are:

```@example
using Entropies, InteractiveUtils
subtypes(ProbabilitiesEstimator)
```

## Probabilities

```@docs
Probabilities
probabilities
probabilities!
ProbabilitiesEstimator
```

## Generalized entropy

```@docs
Entropies.genentropy
```

## Fast histograms

```@docs
Entropies.binhist
```
