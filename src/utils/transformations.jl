


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
