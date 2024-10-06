export DualTotalCorrelation

"""
    DualTotalCorrelation <: MultivariateInformationMeasure
    DualTotalCorrelation(; base = 2)

The dual total correlation measure.

## Compatible estimators

- [`JointProbabilities`](@ref).

## Description 

Let `` \\{X_1, \\ldots, X_n\\} `` be a set of `` n `` random variables.
The dual total correlation `` D(X_1, \\ldots, X_n) `` is given by
```math
D(X_1, \\ldots, X_n) = H(X_1, \\ldots, X_n) - 
\\sum_\\{i=1\\}^\\{n\\} 
H(X_i \\mid X_1, \\ldots, X_\\{i-1\\}, X_\\{i+1\\}, \\ldots, X_n),
```

where `` H(X_1, \\ldots, X_n) `` is the joint entropy of the variables `` \\{X_1, \\ldots, X_n\\} ``,
and `` H(X_i \\mid \\cdots) `` is the conditional entropy of variable `` X_i `` given the rest of 
the variables.
"""
Base.@kwdef struct DualTotalCorrelation{B} <: MultivariateInformationMeasure
    base::B = 2
end

function association(m::DualTotalCorrelation, probs::Probabilities{T, N}) where {T, N}
    mjoint = JointEntropyShannon(; base = m.base)
    mcond = ConditionalEntropyShannon(; base = m.base)
    
end