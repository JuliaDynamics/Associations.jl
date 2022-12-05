# Probabilities

## [Probabilities API](@id probabilities_estimators)

The probabilities API is defined by

- [`ProbabilitiesEstimator`](@ref)
- [`probabilities`](@ref)
- [`probabilities_and_outcomes`](@ref)

```@docs
ProbabilitiesEstimator
probabilities
probabilities!
Probabilities
probabilities_and_outcomes
outcomes
outcome_space
total_outcomes
missing_outcomes
```

## Estimators

### Counting

```@docs
CountOccurrences
```

### Histograms

```@docs
ValueHistogram
RectangularBinning
FixedRectangularBinning
```

### Permutation (symbolic)

```@docs
SymbolicPermutation
SymbolicWeightedPermutation
SymbolicAmplitudeAwarePermutation
SpatialSymbolicPermutation
```

### Dispersion (symbolic)

```@docs
Dispersion
```

### TransferOperator (binning)

```@docs
TransferOperator
```

For explicit estimation of the transfer operator, see
[Entropies.jl](https://github.com/JuliaDynamics/Entropies.jl).

### Kernel density

```@docs
NaiveKernel
```

### Local likelihood

```@docs
LocalLikelihood
```

### Timescales

```@docs
WaveletOverlap
PowerSpectrum
```

### Diversity

```@docs
Diversity
```