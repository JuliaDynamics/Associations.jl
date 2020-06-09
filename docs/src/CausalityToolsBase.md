# CausalityToolsBase.jl

## Discretization 

The following concrete subtypes of `BinningScheme` indicate to algorithms which type of 
partition should be used.

```@docs 
RectangularBinning
```

## Embedding optimization

These convenience types are accepted by some algorithms that use delay embeddings.

```@docs
OptimiseDim
OptimiseDelay
```