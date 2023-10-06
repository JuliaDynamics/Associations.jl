using Test
using CausalityTools
using Random
rng = Xoshiro(1234)

x = StateSpaceSet(rand(rng, 50, 3))
y = StateSpaceSet(rand(rng, 50, 3))
z = StateSpaceSet(rand(rng, 50, 2))
w = rand(rng, ['a', 'b'], 50)
o2 = OrdinalPatternEncoding(2)
o3 = OrdinalPatternEncoding(3)
ow = UniqueElementsEncoding(w)

# Using a single encoding should apply the encoding to all input datasets.
@test probabilities(CodifyPoints(o3), x) isa Probabilities{<:Integer, 1}
@test probabilities(CodifyPoints(o3), x, x) isa Probabilities{<:Integer, 2}

# Using multiple encodings, the number of input encodings must match the number of
# input datasets.
@test probabilities(CodifyPoints(o3, ow), x, w) isa Probabilities{<:Integer, 2}
@test probabilities(CodifyPoints(o3, o3), x, x) isa Probabilities{<:Integer, 2}
@test probabilities(CodifyPoints(o2, o3), z, x) isa Probabilities{<:Integer, 2}
@test probabilities(CodifyPoints(o2, o3, o3), z, x, y) isa Probabilities{<:Integer, 3}

# Length-2 encoding won't work on state vectors of length 3
@test_throws ArgumentError probabilities(CodifyPoints(o2), x)

# When multiple encodings are provided, then the length of the encoding must match
# the length of the points. Here, we accidentally mixed the order of the encodings.
@test_throws ArgumentError probabilities(CodifyPoints(o3, o2, o3), z, x, y)
