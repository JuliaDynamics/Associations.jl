export MIShannon

"""
    MIShannon <: MutualInformation
    MIShannon(; base = 2, definition = MIDefinitionShannonH3())

The Shannon mutual information measure, computed according to `definition` to the
given `base`.

## Description

The Shannon mutual information between two discrete random variables ``X`` and ``Y`` with
supports ``\\mathcal{X}`` and ``\\mathcal{Y}`` is defined as

```math
I(X; Y) = \\sum_{x \\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y)
log \\left( \\dfrac{p(x, y)}{p(x)p(y)} \\right).
```

## Definitions

Estimating ``I(X; Y)`` from data can be done in many different ways, because
there are many possible ways of formulating it. We offer the following
definitions/estimation methods.

- The double-sum definition ([`MIDefinitionShannonDoubleSum`](@ref)), as
    described above.
- The 3-entropies decomposition ([`MIDefinitionShannonH3`](@ref)).

More variants are possible. Pull requests are welcome!

Using a [`MutationalInformationEstimator`](@ref) overrides `definition`, because
the estimator controls the definition and how it is computed.

See also: [`mutualinfo`](@ref).
"""
struct MIShannon{E <: Shannon, D} <: MutualInformation{E, D}
    e::E
    definition::D
    function MIShannon(; base::Real = 2, definition::D = MIDefinitionShannonH3()) where D
        e = Shannon(; base)
        new{typeof(e), D}(e, definition)
    end
    function MIShannon(e::Shannon; definition::D = MIDefinitionShannonH3()) where D
        new{typeof(e), D}(e, definition)
    end
end
