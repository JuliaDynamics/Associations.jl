export TERenyi

"""
    TERenyi() <: TransferEntropy

The Rényi transfer entropy from Jizba et al. (2012).
"""
struct TERenyi{E <: Renyi} <: TransferEntropy
    e::E
    function TERenyi(; base = 2, q = 1.5)
        e = Renyi(; base = base, q = q)
        return new{typeof(e)}(e)
    end
    function TERenyi(e::E) where E <: Renyi
        return new{E}(e)
    end
end

"""
    escort_distribution(probs, i::Int, q::Real)

The escort distribution for a probability distribution `probs`. For `q > 1`, the
escort distribution emphasises more probable events and de-emphasises more improbable
events. For `q < 1`, the situation is reversed.

```math
\\text{esc}_q(x) = \\dfrac{p^q(x)}{\\sum_{x \\in \\mathcal{X}} p^q(x)}
```
"""
function escort_distribution(probs, i::Int, q)
    return probs[i]^q / sum(probs .^ q)
end
