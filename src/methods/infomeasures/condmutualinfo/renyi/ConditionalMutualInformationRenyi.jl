export CMIRenyi

"""
    CMIRenyi <: ConditionalMutualInformation
    CMIRenyi(; base = 2)

`CMIRenyi` a directive to compute the Renyi conditional
mutual information (CMI) between random variables ``X \\in \\mathbb{R}^{d_X}`` and
``Y \\in \\mathbb{R}^{d_Y}``, given ``Z \\in \\mathbb{R}^{d_Z}``.

## Supported definitions

If used with a [`ConditionalMutualInformationEstimator`](@ref), then CMI is estimated
using some direct method, and the docstring for the estimator explains the estimation
procedure and formula.

See also: [`condmutualinfo`](@ref).
"""
struct CMIRenyi{E <: Renyi} <: ConditionalMutualInformation
    e::E
    function CMIRenyi(; base = 2, q = 1.5)
            e = Renyi(; base, q)
        new{typeof(e)}(e)
    end
end
