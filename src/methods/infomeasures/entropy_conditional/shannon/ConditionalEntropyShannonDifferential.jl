export ConditionalEntropyShannonDifferential
export Shannon2h

"""
    Shannon2h <: ConditionalEntropyDefinition

`Shannon2h` is a directive for [`ConditionalEntropyShannonDifferential`](@ref) used
in [`entropy_conditional`](@ref) with an [`EntropyEstimator`](@ref)` to compute the
continuous Shannon conditional entropy using the formula ``h(Y | X) = h(X,Y) - h(X)``.
"""
struct Shannon2h <: ConditionalEntropyDefinition end

"""
    ConditionalEntropyShannonDifferential <: ConditionalEntropyEstimator
    ConditionalEntropyShannonDifferential(; base = 2, definition::Definition = Shannon2h())

`ConditionalEntropyShannon` is a generic plug-in estimator for the continuous conditional
Shannon entropy ``h(Y | X)``.

It computes the discrete conditional entropy to the given `base` by first approximating
probabilities using `est`, and plugging them into the formula given by `definition`.
With default settings, it computes ``h(Y | X) = h(X,Y) - h(X)``

## Supported definitions

- [`Shannon2h`](@ref).

## Usage

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
est = Kraskov(k = 4)
measure = ConditionalEntropyShannonDifferential(; base = 2)
entropy_conditional(measure, est, x, y)
```
See also: [`mutualinfo`](@ref).
"""
struct ConditionalEntropyShannonDifferential{D <: Definition, E <: Renyi} <: ConditionalEntropyEstimator
    e::E
    definition::D
    function ConditionalEntropyShannonDifferential(; base = 2,
            definition::D = Shannon2h()) where {D}
            e = Shannon(; base)
        new{D, typeof(e)}(e, definition)
    end
end

function estimate(measure::ConditionalEntropyShannonDifferential{<:Shannon2h}, est,
        x::AbstractDataset, y::AbstractDataset)
    XY = Dataset(x, y)
    X = Dataset(x)
    return entropy(measure.e, est, XY) - entropy(measure.e, est, X)
end
