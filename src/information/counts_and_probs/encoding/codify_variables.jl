using ComplexityMeasures
import ComplexityMeasures: codify
import ComplexityMeasures: OutcomeSpace

using DelayEmbeddings: embed
export CodifyVariables
export codify

# TODO: implement this Generically for `Encodings` too (will require type-parameterized
# number of elements for the Encodings).

"""
    CodifyVariables <: Discretization
    CodifyVariables(outcome_space::OutcomeSpace)

The `CodifyVariables` discretization scheme quantises input data in a column-wise manner
using the given `outcome_space`.

## Compatible outcome spaces

- [`UniqueElements`](@ref) (for when data are pre-discretized)
- [`BubbleSortSwaps`](@ref)
- [`CosineSimilarityBinning`](@ref)
- [`OrdinalPatterns`](@ref)
- [`Dispersion`](@ref)

# Description

The main difference between `CodifyVariables` and [`CodifyPoints`] is that the former
uses [`OutcomeSpace`](@ref)s for discretization. This usually means that some
transformation is applied to the data before discretizing. For example, some outcome
constructs a delay embedding from the input (and thus encodes sequential information)
before encoding the data.

Specifically, given `x::AbstractStateSpaceSet...`, where the `i`-th dataset `x[i]` 
is assumed to represent a single series of measurements, `CodifyVariables` encodes
 `x[i]` by [`codify`](@ref)-ing into a series of integers 
using an appropriate  [`OutcomeSpace`](@ref). This is typically done by first 
sequentially transforming the data and then running sliding window (the width of 
the window is controlled by `outcome_space`) across the data, and then encoding the 
values within each window to an integer.

## Examples

```julia
using CausalityTools
x, y = rand(100), rand(100)
d = CodifyVariables(OrdinalPatterns(m=2))
cx, cy = codify(d, x, y)
```
"""
struct CodifyVariables{N} <: Discretization{N}
    outcome_spaces::NTuple{N, OutcomeSpace}
    function CodifyVariables(outcome_spaces::NTuple{N, OutcomeSpace}) where N
        if N > 1
            s = "It is currently only possible to use the same `OutcomeSpace` for all " *
                "variables. Got $N different encodings"
            throw(ArgumentError(s))
        end
        new{N}(outcome_spaces)
    end
end

function CodifyVariables(o::OutcomeSpace)
    return CodifyVariables((o,))
end

"""
    codify(d::CodifyVariables, x::Vararg{<:AbstractStateSpaceSet, N})
    codify(d::CodifyPoints, x::Vararg{<:AbstractStateSpaceSet, N})

Codify each timeseries `xᵢ ∈ x` according to the given encoding/discretization `d`.

## Compatible discretizations

- [`CodifyVariables`](@ref)
- [`CodifyPoints`](@ref)

## Examples

```julia
using CausalityTools

# Sliding window encoding
x = [0.1, 0.2, 0.3, 0.2, 0.1, 0.0, 0.5, 0.3, 0.5]
xc1 = codify(CodifyVariables(OrdinalPatterns(m=2)), x) # should give [1, 1, 2, 2, 2, 1, 2, 1]
xc2 = codify(OrdinalPatterns(m=2), x) # equivalent
length(xc1) < length(x) # should be true, because `OrdinalPatterns` delay embeds.  

# Point-by-point encoding
x, y = StateSpaceSet(rand(100, 3)), StateSpaceSet(rand(100, 3))
cx, cy = codify(CodifyPoints(OrdinalPatternEncoding(3)), x, y)
```
"""
function codify(encoding::CodifyVariables, x) end

function codify(encoding::CodifyVariables{1}, x::Vararg{Any, 1})
    e = first(encoding.outcome_spaces)
    x̂ = ComplexityMeasures.codify(e, first(x))
    return x̂::Vector{<:Integer}
end

function codify(encoding::CodifyVariables{1}, x::Vararg{Any, N}) where N
    e = first(encoding.outcome_spaces)
    x̂ = map(xᵢ -> ComplexityMeasures.codify(e, xᵢ), x)
    return x̂::NTuple{N, Vector{<:Integer}}
end

function codify(encoding::CodifyVariables{1}, x::AbstractStateSpaceSet)
    return codify(encoding, columns(x))
end
