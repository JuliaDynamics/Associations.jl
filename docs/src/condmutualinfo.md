# Conditional mutual information

## API

```@docs
condmutualinfo
```

## Types of conditional mutual information

Analogously to [`EntropyDefinition`](@ref), there are also different types of
CMIs (subtypes of [`ConditionalMutualInformation`](@ref)).

```@docs
ConditionalMutualInformation
CMIShannon
CMIRenyi
```

### Definitions

CMIs differ from basic entropies in that there are *multiple* possible
definitions of a CMI. Any of the definitions below can be given
as a keyword to the [`ConditionalMutualInformation`](@ref)s listed above (e.g.
[`CMIShannon`](@ref) or
[`CMIRenyi`](@ref)), indicating how estimation should be done.

```@docs
ConditionalMutualInformationDefinition
CMIDefinitionShannonH4
CMIDefinitionShannonMI2
CMIDefinitionRenyiH4
CMIDefinitionRenySarbu
```

## Estimators

The following estimators compute CMI in some dedicated way, and always
override any of the definitions above.

```@docs
ConditionalMutualInformationEstimator
FrenzelPompeVelmejkaPalus
MesnerShalisi
PoczosSchneiderCMI
Rahimzamani
```
