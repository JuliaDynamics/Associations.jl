export MITsallis

"""
    MITsallis <: MutualInformation
    MITsallis(; base = 2, q = 1.5, definition = MIDefinitionTsallisH3Furuichi())

The [`Tsallis`](@ref) mutual information measure, computed according to `definition` to the
given `base`.

## Definition

Multiple definitions of Tsallis mutual information exist, and their properties and
interpretations differ. We currently implement

- Furuichi (2006)'s discrete definition ([`MIDefinitionTsallisH3Furuichi`](@ref)).
- Martin et al. (2004)'s discrete definition ([`MIDefinitionTsallisH3Martin`](@ref))

More variants of the discrete Tsallis mutual information are possible. We currently don't
implement any continuous variants of the Tsallis mutual information.
Pull requests are welcome!

If used with a [`MutationalInformationEstimator`](@ref), then `definition`
is ignored, and the estimator controls the computed quantity.

## Examples

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(MITsallis(), est, x, y) # default is `MIDefinitionTsallisH3Furuichi` definition

# specifying definition
mutualinfo(MITsallis(definition = MIDefinitionTsallisH3Furuichi()), est, x, y)
mutualinfo(MITsallis(definition = MIDefinitionTsallisH3Martin()), est, x, y)
```

See also: [`mutualinfo`](@ref).
"""
struct MITsallis{E <: Tsallis, D} <: MutualInformation{E, D}
    e::E
    definition::D
    function MITsallis(; definition::D = MIDefinitionTsallisH3Furuichi(), q = 1.5, base = 2) where D
        e = Tsallis(; q, base)
        new{typeof(e), D}(e, definition)
    end
    function MITsallis(e::Tsallis; definition::D = MIDefinitionTsallisH3Furuichi()) where D
        new{typeof(e), D}(e, definition)
    end
end
