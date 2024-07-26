using Test
using CausalityTools
using Random
rng = Xoshiro(1234)

@test_throws ArgumentError CodifyPoints() # need at least one input encoding

c = Counts([2 2; 3 3])
@test marginal(c, dims = 1).cts == [4, 6]
@test marginal(c, dims = 2).cts == [5, 5]
@test marginal(c, dims = 1:2) == c

# No discretization scheme == UniqueElements
@test counts(rand(1:5, 10)) isa Counts{T, 1} where T
# Must have equal lengths
@test_throws ArgumentError c = counts(rand(1:5, 10), rand(1:3, 5))

# Analytic test
x = [1, 2, 3, 1, 2, 3, 1, 2, 3]
y = [1, 2, 1, 2, 1, 2, 1, 2, 1]

@test counts(x, y) == counts(UniqueElements(), x, y)

# ----------------------------
# With `Encodings` directly
# ----------------------------

x = StateSpaceSet(rand(rng, 50, 3))
y = StateSpaceSet(rand(rng, 50, 3))
z = StateSpaceSet(rand(rng, 50, 2))
w = rand(rng, ['a', 'b'], 50)
o2 = OrdinalPatternEncoding(2)
o3 = OrdinalPatternEncoding(3)
ow = UniqueElementsEncoding(w)

# Number of encodings and input datasets must match
d = CodifyPoints(o3, ow)
@test_throws ArgumentError codify(d, x, y, z)

# Using a single encoding should apply the encoding to all input datasets.
@test counts(CodifyPoints(o3), x) isa Counts{<:Integer, 1}
@test counts(CodifyPoints(o3), x, x) isa Counts{<:Integer, 2}

# Using multiple encodings, the number of input encodings must match the number of
# input datasets.
@test counts(CodifyPoints(o3, ow), x, w) isa Counts{<:Integer, 2}
@test counts(CodifyPoints(o3, o3), x, x) isa Counts{<:Integer, 2}
@test counts(CodifyPoints(o2, o3), z, x) isa Counts{<:Integer, 2}
@test counts(CodifyPoints(o2, o3, o3), z, x, y) isa Counts{<:Integer, 3}

# Length-2 encoding won't work on state vectors of length 3
@test_throws ArgumentError counts(CodifyPoints(o2), x)

# When multiple encodings are provided, then the length of the encoding must match
# the length of the points. Here, we accidentally mixed the order of the encodings.
@test_throws ArgumentError counts(CodifyPoints(o3, o2, o3), z, x, y)

x = StateSpaceSet(rand(rng, 10, 2)); 
y = StateSpaceSet(rand(rng, 10, 2));
d_row = CodifyPoints(OrdinalPatternEncoding{2}()); 
@test counts(d_row, x, y) isa Counts{Int, 2}



# Multiple outcome spaces with the same cardinality with `CodifyVariables`
using Test
using CausalityTools
using Random; rng = Xoshiro(1234)
x, y = rand(rng, 100), rand(rng, 100)
# We must use outcome spaces with the same number of total outcomes.
ox = CosineSimilarityBinning(nbins = factorial(3))
oy = OrdinalPatterns(m = 3)

# Now estimate mutual information
discretization = CodifyVariables((ox, oy))
@test counts(discretization, x, y) isa Counts{Int, 2}

# must have same number of outcomes
@test_throws ArgumentError CodifyVariables(OrdinalPatterns(m = 3), OrdinalPatterns(m = 4))