# Transfer entropy

The transfer entropy API is made up of the following functions and types,
which are listed below:

- [`transferentropy`](@ref).
- [`TransferEntropy`](@ref), and its subtypes.
- [`EmbeddingTE`](@ref), which exists to provide embedding instructions to
    subtypes of [`TransferEntropy`](@ref).
- [`TransferEntropyEstimator`](@ref), and its subtypes.

## API

```@docs
transferentropy
EmbeddingTE
optimize_marginals_te
```

## Definitions

```@docs
TransferEntropy
TEShannon
TERenyiJizba
```

## Estimators

```@docs
TransferEntropyEstimator
Zhu1
Lindner
```

## Convenience

### Symbolic transfer entropy

```@docs
SymbolicTransferEntropy
```

### Phase/amplitude transfer entropy

```@docs
Hilbert
Phase
Amplitude
```
