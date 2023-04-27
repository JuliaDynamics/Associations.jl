
"""
    rank_transformation(x::AbstractVector)
    rank_transformation(x::AbstractStateSpaceSet) → ranks::NTuple{D, Vector}

Rank-transform each variable/column of the length-`n` `D`-dimensional StateSpaceSet `x` and return the
rank-transformed variables as a `D`-tuple of length-`n` vectors.

Returns the unscaled `ranks`. Divide by `n` to get an *approximation* to the
empirical cumulative distribution function (ECDF)  `x`.

## Description

Modulo division by `n`, `rank_transformation` does *roughly* the same as naively computing the ECDF as
```julia
[count(xᵢ .<= x)  for xᵢ in x] / length(x)
```

but an order of magnitude faster and with roughly three orders of magnitude less
allocations. The increased efficiency of this function relative to naively computing the
ECDF is
because it uses sorting of the input data to determine ranks,
arbitrarily breaking ties according to the sorting algorithm. Rank ties can therefore
never occur, and equal values are assigned different but close ranks. To preserve
ties, which you might want to do for example when dealing with
categorical or integer-valued data, use (the much slower) [`empcdf`](@ref).
"""
function rank_transformation(x::AbstractStateSpaceSet)
    s = zeros(Int, length(x)) # re-use for each marginal
    [rank_transformation!(s, xⱼ) for xⱼ in columns(x)]
end

function rank_transformation(x::AbstractVector{T}) where T
    N = length(x)
    s = zeros(Int, N)
    return rank_transformation!(s, x)
end

function rank_transformation!(
        s::AbstractVector{Int},
        x::AbstractVector{T}) where T <: Real
    N = length(x)
    r = zeros(N)
    # Break ties arbitrarily by sorting. This means that ties are broken according to the
    # sorting algorithm used, and equal values are assigned different ranks.
    sortperm!(s, x)
    for j in 1:N
        r[s[j]] = j
    end
    return r
end

"""
    empirical_cdf(x::AbstractVector{<:Real}) → x̄::Vector
    empirical_cdf(x::AbstractStateSpaceSet) → x̄::StateSpaceSet

Rank each sample `xᵢ ∈ x`, rescale it the rank to the interval `[0, 1]` and return
the rescaled ranks `x̄`.
"""
function empcdf(x::AbstractVector)
    F̂ = [count(xᵢ .<= x)  for xᵢ in x] / length(x)
end
empcdf(x::AbstractStateSpaceSet{D, T}) where {D, T} =
    NTuple{D, Vector{eltype(1.0)}}(empcdf(xⱼ) for xⱼ in columns(x))

# # An example worked out by hand.
# X = StateSpaceSet([
#     1 8;
#     2 2;
#     3 6;
#     1 5;
#     2 2;
#     3 1;
#     1 8;
#     2 9;
# ])

# using BenchmarkTools
# using Test
# rank_transform(X)
# @btime rank_transform($X);
# # 2.574 ms (23 allocations: 938.06 KiB)
