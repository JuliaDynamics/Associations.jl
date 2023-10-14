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
encode
decode
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
