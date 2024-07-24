using Test
using Random
rng = Xoshiro(1234)
x, y, z = rand(100), rand(100), rand(100)
@test optimize_marginals_te(OptimiseTraditional(), x, y) isa EmbeddingTE
@test optimize_marginals_te(OptimiseTraditional(), x, y, z) isa EmbeddingTE