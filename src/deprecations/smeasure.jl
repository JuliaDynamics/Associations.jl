
export s_measure

"""
    s_measure(measure::SMeasure, x::VecOrSSSet, y::VecOrSSSet) → s ∈ [0, 1]

Compute the given `measure` to quantify the directional dependence between
univariate/multivariate time series `x` and `y`.

Returns a scalar `s` where `s = 0` indicates independence between `x` and `y`,
and higher values indicate synchronization between `x` and `y`, with complete
synchronization for `s = 1.0`.

## Example

```julia
using CausalityTools

# A two-dimensional Ulam lattice map
sys = ulam(2)

# Sample 1000 points after discarding 5000 transients
orbit = trajectory(sys, 1000, Ttr = 5000)
x, y = orbit[:, 1], orbit[:, 2]

# 4-dimensional embedding for `x`, 5-dimensional embedding for `y`
m = SMeasure(dx = 4, τx = 3, dy = 5, τy = 1)
s_measure(m, x, y)
```

"""
function s_measure(x::VecOrSSSet, y::VecOrSSSet; kwargs...)
    if !isempty(kwargs)
        @warn(
            "Providing keywords to `s_measure` is deprecated. " *
            "Use `s_measure(SMeasure(; kwargs...), source, target) instead of " *
            "`s_measure(source, target; kwargs...)`"
        )
    end
    return estimate(SMeasure(; kwargs...), x, y)
end
