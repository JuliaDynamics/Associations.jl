
# Probabilities API

The probabilities API is defined by

- [`ProbabilitiesEstimator`](@ref), and its subtypes.
- [`probabilities`](@ref)
- [`probabilities_and_outcomes`](@ref)

See also [contingency tables](@ref contingency_table_api) for a multivariate version.

The probabilities API is re-exported from [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl). Why? Most discrete information theoretic association measures are estimated
using some sort of [`ProbabilitiesEstimator`](@ref)s, because their formulas are simply functions
of probability mass functions.

## Probabilities

```@docs
ProbabilitiesEstimator
probabilities
probabilities!
Probabilities
```

## Estimators

### Overview

Here, we list probabilities estimators that are compatible with CausalityTools.jl. Note that not
all probabilities estimators from ComplexityMeasures.jl are included. This is because for
the information-based association measures here, the probabilities estimator must be
compatible with multivariate data, or have an implementation for [`marginal_encodings`](@ref),
which discretizes each dimension of the multivariate input data separately.

| Estimator                     | Principle                                      |
| :---------------------------- | :--------------------------------------------- |
| [`Contingency`](@ref)         | Count co-occurrences, optionally discretize first |
| [`CountOccurrences`](@ref)    | Count of unique elements                       |
| [`ValueHistogram`](@ref)      | Binning (histogram)                            |
| [`TransferOperator`](@ref)    | Binning (transfer operator)                    |
| [`NaiveKernel`](@ref)         | Kernel density estimation                      |
| [`SymbolicPermutation`](@ref) | Ordinal patterns                               |
| [`Dispersion`](@ref)          | Dispersion patterns                            |

### Contingency

```@docs
Contingency
```

### Count occurrences

```@docs
CountOccurrences
```

### Histograms (binning)

```@docs
ValueHistogram
RectangularBinning
FixedRectangularBinning
```

### Transfer operator (binning)

```@docs
TransferOperator
```

#### Utility methods/types

For explicit estimation of the transfer operator, see
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl).

```@docs
InvariantMeasure
invariantmeasure
transfermatrix
```

### Symbolic permutations

```@docs
SymbolicPermutation
```

### Dispersion patterns

```@docs
Dispersion
```

### Kernel density

```@docs
NaiveKernel
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

## Utilities

### Outcomes

```@docs
probabilities_and_outcomes
outcomes
outcome_space
total_outcomes
missing_outcomes
```

### Encodings

Some probability estimators first "encode" input data into an intermediate representation indexed by the positive integers. This intermediate representation is called an "encoding".

The encodings API is defined by:

- [`Encoding`](@ref)
- [`encode`](@ref)
- [`decode`](@ref)

```@docs
Encoding
encode
decode
```

#### Available encodings

```@docs
OrdinalPatternEncoding
GaussianCDFEncoding
RectangularBinEncoding
```
