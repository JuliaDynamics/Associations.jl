export MIShannon

"""
    MIShannon <: MutualInformation
    MIShannon(; base = 2)

The Shannon mutual information measure.

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
struct MIShannon{E <: Renyi} <: MutualInformation
    e::E
    function MIShannon(; base = 2) where {D}
        e = Shannon(; base)
        new{typeof(e)}(e)
    end
    function MIShannon(e::E) where {E <: Renyi}
        @assert e.q ≈ 1 || error("MIShannon not defined with q = $(e.q)")
        new{E}(e)
    end
end

# Define default behaviour.
estimate(measure::MIShannon, est, x, y) = estimate(ShannonH3(), measure, est, x, y)
