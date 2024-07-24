using Test
using Random
rng = Xoshiro(1234)

# ----------------------------------------------------------------
# Embedding optimizations
# ----------------------------------------------------------------
x, y, z = rand(100), rand(100), rand(100)
@test OptimizeTraditional() isa OptimiseTraditional() # alias
@test optimize_marginals_te(OptimiseTraditional(), x, y) isa EmbeddingTE
@test optimize_marginals_te(OptimiseTraditional(), x, y, z) isa EmbeddingTE
@test optimize_marginals_te(OptimiseTraditional(), x, y, z; exclude_source = true) isa EmbeddingTE
@test EmbeddingTE(OptimiseTraditional(), x, y, z) isa EmbeddingTE
