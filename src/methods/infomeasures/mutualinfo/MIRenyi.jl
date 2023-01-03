export MIRenyi

"""
    MIRenyi <: MutualInformation
    MIRenyi(; base = 2, q = 1.5, definition = MIDefinitionRenyiSarbu())

The [`Renyi`](@ref) mutual information measure, computed according to `definition` to the
given `base`.

## Definition

Multiple definitions of Rényi mutual information exist, and their properties and
interpretations differ. We currently implement

- [`MIDefinitionRenyiSarbu`](@ref)). Discrete Rényi mutual information based on the
    Rényi ``\\alpha``-divergences.

More variants of the discrete Renyi mutual information are possible. We currently don't
implement any continuous variants of the Renyi mutual information.
Pull requests are welcome!

If used with a [`MutationalInformationEstimator`](@ref), then `definition`
is ignored, and the estimator controls the computed quantity.

## Examples

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(MIRenyi(definition = MIDefinitionRenyiSarbu()), est, x, y)
```

See also: [`mutualinfo`](@ref).
"""
struct MIRenyi{E <: Renyi, D} <: MutualInformation{E, D}
    e::E
    definition::D
    function MIRenyi(; definition::D = MIDefinitionRenyiH3Furuichi(), q = 1.5, base = 2) where D
        e = Renyi(; q, base)
        new{typeof(e), D}(e, definition)
    end
    function MIRenyi(e::Renyi; definition::D = MIDefinitionRenyiH3Furuichi()) where D
        new{typeof(e), D}(e, definition)
    end
end
