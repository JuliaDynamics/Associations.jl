# [Probabilities](@id probabilities_estimators)

## Probabilities API

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

## Count occurrences (counting)

```@docs
CountOccurrences
```

## Visitation frequency (histograms)

```@docs
ValueHistogram
RectangularBinning
FixedRectangularBinning
```

## Permutation (symbolic)

```@docs
SymbolicPermutation
SymbolicWeightedPermutation
SymbolicAmplitudeAwarePermutation
SpatialSymbolicPermutation
```

## Dispersion (symbolic)

```@docs
Dispersion
```

## Transfer operator (binning)

```@docs
TransferOperator
```

### Utility methods/types

```@docs
InvariantMeasure
invariantmeasure
transfermatrix
```

## Kernel density

```@docs
NaiveKernel
```

## Local likelihood

```@docs
LocalLikelihood
```

## Timescales

```@docs
WaveletOverlap
PowerSpectrum
```

## Diversity

```@docs
Diversity
```
