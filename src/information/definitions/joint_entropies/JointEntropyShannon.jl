using ComplexityMeasures: Shannon

export JointEntropyShannon

"""
    JointEntropyShannon <: JointEntropy
    JointEntropyShannon(; base = 2)

The Shannon joint entropy measure [CoverThomas1999](@cite).

## Usage 

- Use with [`association`](@ref) to compute the Shannon joint entropy between 
    two variables.

## Compatible estimators

- [`JointProbabilities`](@ref)

## Definition

Given two two discrete random variables ``X`` and ``Y`` with ranges ``\\mathcal{X}`` and
``\\mathcal{X}``, [CoverThomas1999](@citet) defines the Shannon joint entropy as

```math
H^S(X, Y) = -\\sum_{x\\in \\mathcal{X}, y \\in \\mathcal{Y}} p(x, y) \\log p(x, y),
```

where we define ``log(p(x, y)) := 0`` if ``p(x, y) = 0``.

## Examples

```julia
using CausalityTools
x, y = rand(100), rand(100)
measure = JointEntropyShannon()
discretization = CodifyVariables(Dispersion(m=2, c = 3))
est = JointProbabilities(measure, discretization)
information(est, x, y)
```
"""
Base.@kwdef struct JointEntropyShannon{B} <: JointEntropy
    base::B = 2
end

function information(definition::JointEntropyShannon, pxy::Probabilities{T, 2}) where T
    (; base) = definition
    
    h = 0.0
    for p in pxy
        if p != 0 # Define log(0) = 0
            h += p * log(p)
        end
    end
    h = -h
    return _convert_logunit(h, ℯ, base)
end