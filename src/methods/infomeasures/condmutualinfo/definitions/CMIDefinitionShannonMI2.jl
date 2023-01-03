export CMIDefinitionShannonMI2

"""
    CMIDefinitionShannonMI2 <: ConditionalMutualInformationDefinition
    CMIDefinitionShannonMI2()

A directive to compute conditional mutual information (CMI) as a difference of mutual
information terms.

```math
I(X; Y | Z) = I^S(X; Y, Z) + I^S(X; Y)
```

where ``H^S(\\cdot)`` is the Shannon entropy of the argument. The formula applies both
for discrete and continuous CMI.
```
"""
Base.@kwdef struct CMIDefinitionShannonMI2{M <: MutualInformation} <: CMIMI2{M}
    measure::M = MIShannon()
end
