# TransferEntropy.jl 

CausalityTools provides the following interface for computing transfer entropy. The `transferentropy` method 
dispatches on different `TransferEntropyEstimator`s (listed below).

```@docs
transferentropy
```

## [Estimators](@id te_estimators)

```@docs
VisitationFrequency
TransferOperatorGrid
NearestNeighborMI
SimplexEstimator
SymbolicPerm
SymbolicAmplitudeAware
```

## Generalized embedding

Information about the dimensions and lags of the marginals used for 
transfer entropy computation are given by an `EmbeddingTE` instance.

```@docs
EmbeddingTE
```

## Binning heuristics

The binning-based estimators rely on a coarse-graining of the reconstructed 
state space to compute relevant marginal entropies. The folllowing heuristic
methods compute a suitable coarse-graining based on the number of available 
points and the dimension of the reconstructed phase space. 

```@docs
PalusLimit
ExtendedPalusLimit
```