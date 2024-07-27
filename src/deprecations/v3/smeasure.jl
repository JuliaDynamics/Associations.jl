
export s_measure

"""
    s_measure(measure::SMeasure, x::VectorOrStateSpaceSet, y::VectorOrStateSpaceSet) → s ∈ [0, 1]

Compute the given `measure` to quantify the directional dependence between
univariate/multivariate time series `x` and `y`.

Returns a scalar `s` where `s = 0` indicates independence between `x` and `y`,
and higher values indicate synchronization between `x` and `y`, with complete
synchronization for `s = 1.0`.

## Example

```julia
using Associations

x, y = rand(1000), rand(1000)

# 4-dimensional embedding for `x`, 5-dimensional embedding for `y`
m = SMeasure(dx = 4, τx = 3, dy = 5, τy = 1)
association(m, x, y)
```
"""
function s_measure(x, y; kwargs...)
    @warn(
        "Convenience function `s_measure` is deprecated. " *
        "Use `association(SMeasure(; kwargs...), x, y)` instead."
    )
    return association(SMeasure(; kwargs...), x, y)
end

function s_measure(measure::SMeasure, x, y; kwargs...)
    @warn(
        "Convenience function `s_measure` is deprecated. " *
        "Use `association(SMeasure(; kwargs...), x, y)` instead."
    )
    return association(SMeasure(; kwargs...), x, y)
end