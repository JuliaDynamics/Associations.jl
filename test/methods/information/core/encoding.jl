using Test
using CausalityTools
using Random
rng = Xoshiro(1234)

# Constructors
@test_throws ArgumentError CodifyVariables((OrdinalPatterns(m=2), OrdinalPatterns(m=3)))

x = StateSpaceSet(rand(rng, 50, 3))
y = StateSpaceSet(rand(rng, 50, 3))
z = StateSpaceSet(rand(rng, 50, 2))
w = rand(rng, ['a', 'b'], 50)
o2 = OrdinalPatternEncoding(2)
o3 = OrdinalPatternEncoding(3)
ow = UniqueElementsEncoding(w)

# Using a single encoding should apply the encoding to all input datasets.
@test codify(CodifyPoints(o3), x) isa Vector{<:Integer}
@test codify(CodifyPoints(o3), x, x) isa NTuple{2, Vector{<:Integer}}

# Using multiple encodings, the number of input encodings must match the number of
# input datasets.
@test codify(CodifyPoints(o3, ow), x, w) isa NTuple{2, Vector{<:Integer}}
@test codify(CodifyPoints(o3, o3), x, x) isa NTuple{2, Vector{<:Integer}}
@test codify(CodifyPoints(o2, o3), z, x) isa NTuple{2, Vector{<:Integer}}
@test codify(CodifyPoints(o2, o3, o3), z, x, y) isa NTuple{3, Vector{<:Integer}}

# Length-2 encoding won't work on state vectors of length 3
@test_throws ArgumentError codify(CodifyPoints(o2), x)

# When multiple encodings are provided, then the length of the encoding must match
# the length of the points. Here, we accidentally mixed the order of the encodings.
@test_throws ArgumentError codify(CodifyPoints(o3, o2, o3), z, x, y)

#----------------------------------------------------------------
# Per variable/column encoding
#----------------------------------------------------------------

# Single variables
x = rand(rng, 100)
o = ValueBinning(3)
@test codify(CodifyVariables(o), x) isa Vector{<:Integer}
@test codify(CodifyVariables(o), (x, )) isa NTuple{1, Vector{<:Integer}}

# Point-by-point encoding
x, y = StateSpaceSet(rand(100, 3)), StateSpaceSet(rand(100, 3))
cx, cy = codify(CodifyPoints(OrdinalPatternEncoding(3)), x, y)
@test cx isa Vector{Int}
@test cy isa Vector{Int}
