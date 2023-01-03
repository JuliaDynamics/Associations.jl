export CMIRenyi

"""
    CMIRenyi <: ConditionalMutualInformation
    CMIRenyi(; base = 2, definition = CMIDefinitionRenyiSarbu())

`CMIRenyi` a directive to compute the Renyi conditional
mutual information (CMI) between random variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}``, given ``Z \\in \\mathbb{R}^{d_Z}``.

## Supported definitions

- [`CMIDefinitionRenyiSarbu`](@ref). A definition of discrete Rényi CMI based on
    conditional Rényi ``\\alpha``-divergences.

If used with a [`ConditionalMutualInformationEstimator`](@ref), then CMI is estimated
using some direct method, and the docstring for the estimator explains the estimation
procedure and formula.

See also: [`condmutualinfo`](@ref).
"""
struct CMIRenyi{E <: Renyi, D} <: ConditionalMutualInformation{E, D}
    e::E
    definition::D
    function CMIRenyi(; base = 2, q = 1.5, definition::D = CMIDefinitionRenyiSarbu()) where D
            e = Renyi(; base, q)
        new{typeof(e), D}(e, definition)
    end
    function CMIShannon(e::E; definition::D = CMIDefinitionRenyiH4()) where {E <: Shannon, D}
        new{E, D}(e, definition)
    end
end
