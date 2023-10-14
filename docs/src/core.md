# Core


## Discretization

```@docs
Discretization
codify
```

### Column-wise discretization

```@docs
CodifyVariables
```

### Point-wise discretization

The [`CodifyVariables`](@ref) discretization scheme encodes input data points directly 
to integers, *without* applying any transformation to the data before encoding.

```@docs
CodifyPoints
Encoding
encoding
GaussianCDFEncoding
OrdinalPatternEncoding
RelativeMeanEncoding
RelativeFirstDifferenceEncoding
UniqueElementsEncoding
CombinationEncoding
RectangularBinEncoding
```

## Counting and probabilities

For counting and probabilities, CausalityTools.jl extends the single-variable machinery
in ComplexityMeasures.jl to multiple variables.

```@docs
CausalityTools.Counts
CausalityTools.counts
CausalityTools.Probabilities
CausalityTools.probabilities
marginal
```
