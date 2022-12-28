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
Base.@kwdef struct CMIShannon{E <: Shannon} <: ConditionalMutualInformation
    e::E = Shannon()
    function CMIShannon(; base::T = 2) where T <: Real
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
end

# Different ways of estimating the CMI. More are possible.
include("CMI2MIShannon.jl")
include("CMI4HShannon.jl")
