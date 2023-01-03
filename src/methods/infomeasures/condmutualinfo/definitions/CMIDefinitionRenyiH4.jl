export CMIDefinitionRenyiH4

"""
    CMIDefinitionRenyiH4 <: ConditionalMutualInformationDefinition
    CMIDefinitionRenyiH4()

A directive to compute Renyi conditional mutual information (CMI) as a sum of four
entropy terms

```math
I(X; Y | Z) = H^R(X, Z) + H^R(Y, z) - H^R(X, Y, Z) - H^R(Z),
```

where ``H^R(\\cdot)`` is the [`Renyi`](@ref) entropy of the argument. The formula applies for
the continuous case too, just changing notation to lower-case ``h`` to indicate that the
*differential* entropy is used:

```math
I(X; Y | Z) = h^R(X, Z) + h^R(Y, z) - h^R(X, Y, Z) - h^R(Z).
```

!!! info
    Rényi CMO do not share all properties of Shannon CMI.
"""
struct CMIDefinitionRenyiH4 <: CMIH4 end
