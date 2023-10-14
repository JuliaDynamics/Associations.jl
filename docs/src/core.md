# Core

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

## Discretization

There are many ways of discretizing multiple input datasets. The [`CodifyPoints`](@ref) 
discretization scheme, on the other hand, encodes on each `N`-dimensional point
in each input dataset directly as an integer, using some 

```@docs
Discretization
```

### Column-wise discretization

The [`CodifyVariables`](@ref)
discretization scheme quantises input data in a column-wise manner, often applying 
some transformation to the data first, e.g such delay embeddings, and thus 
uses sequential information to encode data. This discretization scheme integrates 
nicely with[`OutcomeSpace`](@ref)s from ComplexityMeasures.jl, which can be used 
interchangeable with [`CodifyVariables`](@ref) for the simplest use cases.

```@docs
CodifyVariables
```

### Point-wise discretization

The [`CodifyVariables`](@ref) discretization scheme encodes input data points directly 
to integers, *without* applying any transformation to the data before encoding.

```@docs
CodifyPoints
ComplexityMeasures.Encoding
ComplexityMeasures.encoding
```