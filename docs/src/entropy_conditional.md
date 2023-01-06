# Conditional entropy

The function that estimates conditional entropy from data is [`entropy_conditional`](@ref).
It can estimate different types of conditional entropies from the scientific literature,
and each is represented here as a subtype of [`ConditionalEntropy`](@ref). Because a
conditional entropy can be formulated in many different ways, each conditional entropy
type can be estimated according to multiple definitiosn, which are represented by
subtypes of [`ConditionalEntropyDefinition`](@ref)).

To see which estimators are compatible with the various definitions, see the
[overview table](@ref conditionalentropy_overview) below.

## Conditional entropy API


```@docs
entropy_conditional
ConditionalEntropy
ConditionalEntropyDefinition
```

## Conditional entropy definitions

More variants of the discrete Tsallis conditional entropy are possible beyond those
defined below.
Pull requests are welcome!

### Shannon conditional entropy

```@docs
CEShannon
```

### Tsallis conditional entropy

```@docs
CETsallisFuruichi
CETsallisAbe
```

## Discrete conditional entropy

### [Table of discrete conditional entropy estimators] (@id conditionalentropy_overview)

Here, we list the [`ProbabilitiesEstimator`](@ref)s that are compatible with
[`entropy_conditional`](@ref), and which definitions they are valid for.

| Estimator                    | Principle           |     Input data     | [`CEShannon`](@ref) | [`CETsallisFuruichi`](@ref) | [`CETsallisAbe`](@ref) |
| ---------------------------- | ------------------- | :----------------: | :-----------------: | :-------------------------: | :--------------------: |
| [`CountOccurrences`](@ref)   | Frequencies         | `Vector`/`Dataset` |         ✅          |             ✅              |           ✅           |
| [`ValueHistogram`](@ref)     | Binning (histogram) | `Vector`/`Dataset` |         ✅          |             ✅              |           ✅           |
| [`SymbolicPermuation`](@ref) | Ordinal patterns    |      `Vector`      |         ✅          |             ✅              |           ✅           |
| [`Dispersion`](@ref)         | Dispersion patterns |      `Vector`      |         ✅          |             ✅              |           ✅           |

## Differential conditional entropy

We currently don't support the computation of differential conditional entropy.
Pull requests are welcome!
