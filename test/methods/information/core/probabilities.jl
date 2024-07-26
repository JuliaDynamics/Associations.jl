using Test
using CausalityTools
using Random
rng = Xoshiro(1234)

x = StateSpaceSet(rand(rng, 50, 3))
y = StateSpaceSet(rand(rng, 50, 3))
z = StateSpaceSet(rand(rng, 50, 2))
w = rand(rng, ['a', 'b'], 50)
o2 = OrdinalPatternEncoding{2}()
o3 = OrdinalPatternEncoding{3}()
ow = UniqueElementsEncoding(w)

# Using a single encoding should apply the encoding to all input datasets.
@test probabilities(CodifyPoints(o3), x) isa Probabilities{T, 1} where T
@test probabilities(CodifyPoints(o3), x, x) isa Probabilities{T, 2} where T

# Using multiple encodings, the number of input encodings must match the number of
# input datasets.
@test probabilities(CodifyPoints(o3, ow), x, w) isa Probabilities{T, 2} where T
@test probabilities(CodifyPoints(o3, o3), x, x) isa Probabilities{T, 2} where T
@test probabilities(CodifyPoints(o2, o3), z, x) isa Probabilities{T, 2} where T
@test probabilities(CodifyPoints(o2, o3, o3), z, x, y) isa Probabilities{T, 3} where T

# Length-2 encoding won't work on state vectors of length 3
@test_throws ArgumentError probabilities(CodifyPoints(o2), x)

# When multiple encodings are provided, then the length of the encoding must match
# the length of the points. Here, we accidentally mixed the order of the encodings.
@test_throws ArgumentError probabilities(CodifyPoints(o3, o2, o3), z, x, y)


# ----------------------------------------------------------------
# Multiple outcome spaces with CodifyVariables
# ----------------------------------------------------------------
x, y = rand(rng, 100), rand(rng, 100)
# We must use outcome spaces with the same number of total outcomes.
ox = CosineSimilarityBinning(nbins = factorial(3))
oy = OrdinalPatterns(m = 3)

# Now estimate mutual information
discretization = CodifyVariables((ox, oy))
@test probabilities(discretization, x, y) isa Probabilities{T, 2} where T