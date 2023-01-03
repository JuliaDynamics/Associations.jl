export CMIShannon

"""
    CMIShannon <: ConditionalMutualInformation
    CMIShannon(; base = 2)

`CMIShannon` a directive to compute the Shannon conditional
mutual information (CMI) between random variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}``, given ``Z \\in \\mathbb{R}^{d_Z}``.

## Supported definitions

If used with a [`ConditionalMutualInformationEstimator`](@ref), then CMI is estimated
using some direct method, and the docstring for the estimator explains the estimation
procedure and formula.

We also support other definitions based on lower-level information measures:

- [`CMI4HShannon`](@ref)
- [`CMI2MIShannon`](@ref)

See also: [`condmutualinfo`](@ref).
"""
struct CMIShannon{E <: Shannon, D} <: ConditionalMutualInformation{E, D}
    e::E
    definition::D
    function CMIShannon(; base::T = 2, definition::D = CMIDefinitionShannonH4()) where {T <: Real, D}
        e = Shannon(; base)
        new{typeof(e), D}(e, definition)
    end
    function CMIShannon(e::E; definition::D = CMIDefinitionShannonH4()) where {E <: Shannon, D}
        new{E, D}(e, definition)
    end
end

# This is necessary because if `MutualInformationEstimator` are used, then they
# require the `measure.definition` keyword, which doesn't exist for all definitions.
function default_measure(measure::CMI{<:Shannon, D}, est::MutualInformationEstimator) where D
    CMIShannon(measure.e, definition = CMIDefinitionShannonMI2())
end
