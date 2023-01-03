export CMIDefinitionShannonH4

"""
    CMIDefinitionShannonH4 <: ConditionalMutualInformationDefinition
    CMIDefinitionShannonH4()

A directive to compute Shannon conditional mutual information (CMI) as a sum of four
entropy terms

```math
I(X; Y | Z) = H^S(X, Z) + H^S(Y, z) - H^S(X, Y, Z) - H^S(Z),
```

where ``H^S(\\cdot)`` is the generalized entropy of the argument. The formula applies for
the continuous case too, just changing notation to lower-case ``h`` to indicate that the
*differential* entropy is used:

```math
I(X; Y | Z) = h^S(X, Z) + h^S(Y, z) - h^S(X, Y, Z) - h^G(Z).
```
"""
struct CMIDefinitionShannonH4 <: CMIH4 end
