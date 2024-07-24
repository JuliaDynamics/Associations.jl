import ComplexityMeasures: codify
import ComplexityMeasures: counts
using ComplexityMeasures: Encoding

export CodifyPoints
export codify

"""
    CodifyPoints{N}
    CodifyPoints(encodings::NTuple{N, Encoding})

`CodifyPoints` points is a [`Discretization`](@ref) scheme that encodes input data points
*without* applying any sequential transformation to the input (as opposed to 
[`CodifyVariables`](@ref), which may apply some transformation before encoding).

## Usage

- Use with [`codify`](@ref)` to encode/discretize input variable on a point-by-point basis.

## Compatible encodings

- [`GaussianCDFEncoding`](@ref)
- [`OrdinalPatternEncoding`](@ref)
- [`RelativeMeanEncoding`](@ref)
- [`RelativeFirstDifferenceEncoding`](@ref)
- [`UniqueElementsEncoding`](@ref)
- [`RectangularBinEncoding`](@ref)
- [`CombinationEncoding`](@ref)

## Description

Given `x::AbstractStateSpaceSet...`, where the `i`-th dataset is assumed to represent
a single series of measurements, `CodifyPoints` encodes each point `pₖ ∈ x[i]` 
using some [`Encoding`](@ref)(s), *without* applying any (sequential) transformation to
the `x[i]` first. This behaviour is different to [`CodifyVariables`](@ref), which
*does* apply a transformation to `x[i]` before encoding.

If `length(x) == N` (i.e. there are `N` input dataset), then `encodings` must be a tuple
of `N` [`Encoding`](@ref). Alternatively, if `encodings` is a single [`Encoding`](@ref),
then that same encoding is applied to every `x[i]`.

## Examples

```julia
using CausalityTools

# The same encoding on two input datasets
x = StateSpaceSet(rand(100, 3))
y = StateSpaceSet(rand(100, 3))
encoding_ord = OrdinalPatternEncoding(3)
cx, cy = codify(CodifyPoints(encoding_ord), x, y)

# Different encodings on multiple datasets
z = StateSpaceSet(rand(100, 2))
encoding_bin = RectangularBinEncoding(RectangularBinning(3), z)
d = CodifyPoints(encoding_ord, encoding_ord, encoding_bin)
cx, cy, cz = codify(d, x, y, z)
```
"""
struct CodifyPoints{N} <: Discretization{N}
    encodings::NTuple{N, Encoding}
    function CodifyPoints(encodings::NTuple{N, Encoding}) where N
        if !(N ≥ 1)
            throw(ArgumentError("CodifyPoints requires at least 1 dimensions"))
        end
        new{N}(encodings)
    end
end
Base.getindex(e::CodifyPoints, i) = getindex(e.encodings, i)

function CodifyPoints(encodings::Vararg{Encoding, N}) where N
    return CodifyPoints(tuple(encodings...))
end

"""
    codify(encoding::CodifyPoints{N}, x::Vararg{<:AbstractStateSpaceSet, N})

Codify each timeseries `xᵢ ∈ x` according to the given `encoding`.

## Examples

```julia
x = StateSpaceSet(rand(10000, 2))
y = StateSpaceSet(rand(10000, 3))
z = StateSpaceSet(rand(10000, 2))

# For `x`, we use a relative mean encoding.
ex = RelativeMeanEncoding(0.0, 1.0, n = 3)
# For `y`, we use a combination encoding.
ey = CombinationEncoding(
    RelativeMeanEncoding(0.0, 1.0, n = 2), 
    OrdinalPatternEncoding(3)
)
# For `z`, we use ordinal patterns to encode.
ez = OrdinalPatternEncoding(2)

# Codify two input datasets gives a 2-tuple of Vector{Int}
codify(CodifyPoints(ex, ey), x, y)

# Codify three input datasets gives a 3-tuple of Vector{Int}
codify(CodifyPoints(ex, ey, ez), x, y, z)
```
"""
function codify(encoding::CodifyPoints, x) end

function codify(encoding::CodifyPoints{1}, x::Vararg{Any, 1})
    e = first(encoding.encodings)
    x̂ = codify_individual_dataset(e, first(x))
    return x̂::Vector{<:Integer}
end

# Apply the same encoding to all input datasets.
function codify(encoding::CodifyPoints{1}, x::Vararg{Any, M}) where {M}
    verify_input(encoding, x...)
    e = first(encoding.encodings)
    x̂ = map(k -> codify_individual_dataset(e, x[k]), tuple(1:M...))

    return x̂::NTuple{M, Vector{<:Integer}}
end


function codify(encoding::CodifyPoints{N}, x::Vararg{Any, M}) where {N, M}
    verify_input(encoding, x...)
    x̂ = map(k -> codify_individual_dataset(encoding[k], x[k]), tuple(1:M...))

    return x̂::NTuple{M, Vector{<:Integer}}
end

function verify_input(encoding::CodifyPoints{N}, x...) where N
    M = length(x)
    if N != M && N != 1
        s = "The given `encoding` is for $N input datasets. $M input datasets were given."
        throw(ArgumentError(s))
    end
    Ls = length.(x)
    if !allequal(Ls)
        throw(ArgumentError("All input datasets must have the same length."))
    end
end

function codify_individual_dataset(encoding::Encoding, x)
    if !(typeof(x) <: AbstractStateSpaceSet)
        encoding = UniqueElementsEncoding(x)
        x̂ = encode.(Ref(encoding), x)
        return x̂
    end

    # x̂[i] := the integer code for the state vector `x[i]`.
    x̂ = zeros(Int, length(x))
    @inbounds for i in eachindex(x)
        x̂[i] = encode(encoding, x[i])
    end
    return x̂
end

 # The decoding step on the second-to-last line is not possible without actually providing
 # the encodings. Therefore, we need to override the Generic implementation of
 # `counts`.
function counts(encoding::CodifyPoints, x...)
    # This converts each dataset `x[i]::StateSpaceSet` into `x̂[i]::Vector{Int}`,
    # where `length(x[i]) == length(x̂[i])`.
    x̂ = codify(encoding, x...)
    # lmaps[i]: a `Dict{outcome_type, Int}` containing the conversion between the
    #   internally encoded outcomes for the `i`-th input, and the actual outcomes
    #   for the `i`-th input.
    cts, lmaps, encoded_outcomes = counts_table(x̂...)

    # Actual outcomes (these outcomes map correspond to those in `x̂`).
    # We can't actually decode any further than this.
    L = length(x)
    outcomes = map(i -> to_outcomes(lmaps[i], encoded_outcomes[i]), tuple(1:L...))

    # Marginal labels are the decoded outcomes.
    decoded_outcomes = map(i -> decode_outcomes(encoding[i], outcomes[i]), tuple(1:L...))
    return Counts(cts, decoded_outcomes)
end

function decode_outcomes(encoding::Encoding, outcomes::Vector{<:Integer})
    return ComplexityMeasures.decode.(Ref(encoding), outcomes)
end
