
# Probabilities API

The probabilities API is defined by

- [`ProbabilitiesEstimator`](@ref), and its subtypes.
- [`probabilities`](@ref)
- [`probabilities_and_outcomes`](@ref)

See also [contingency tables](@ref contingency_table_api) for a multivariate version.

The probabilities API is re-exported from [ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl). Why? Most discrete information theoretic association measures are estimated
using some sort of [`ProbabilitiesEstimator`](@ref)s, because their formulas are simply functions
of probability mass functions.

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

### Kernel density

```@docs
NaiveKernel
```

### Timescales

```@docs
WaveletOverlap
PowerSpectrum
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
