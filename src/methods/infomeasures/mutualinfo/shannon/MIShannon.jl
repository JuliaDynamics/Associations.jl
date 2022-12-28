export MIShannon

"""
    MIShannon <: MutualInformation
    MIShannon(; base = 2)

The Shannon mutual information measure. Used with [`mutualinfo`](@ref).

## Supported definitions

- [`ShannonH3`](@ref).

If used with a [`MutationalInformationEstimator`](@ref), then the estimator dictates the
definition (i.e. which formula is computed).

## Usage

```julia
using CausalityTools
x, y = Dataset(rand(1000)), Dataset(rand(1000))
b = FixedRectangularBinning(0, 1, 5)
estimate(MIShannon(), ValueHistogram(b), x, y)
```
See also: [`mutualinfo`](@ref).
"""
Base.@kwdef struct MIShannon{E <: Shannon} <: MutualInformation
    e::E = Shannon()
    function MIShannon(; base::T = 2) where T <: Real
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
end

# Define default behaviour.
estimate(measure::MIShannon, est, x, y) = estimate(ShannonH3(), measure, est, x, y)

include("ShannonH3.jl")
