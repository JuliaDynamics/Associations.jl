using Test
using Random
rng = Xoshiro(1234)

# ----------------------------------------------------------------
# Embedding optimizations
# ----------------------------------------------------------------
x, y, z = rand(100), rand(100), rand(100)
@test typeof(OptimizeTraditional()) == typeof(OptimiseTraditional()) # alias
@test optimize_marginals_te(OptimiseTraditional(), x, y; exclude_source = true) isa EmbeddingTE
@test optimize_marginals_te(OptimiseTraditional(), x, y; exclude_source = false) isa EmbeddingTE
@test optimize_marginals_te(OptimiseTraditional(), x, y, z; exclude_source = true) isa EmbeddingTE
@test optimize_marginals_te(OptimiseTraditional(), x, y, z; exclude_source = false) isa EmbeddingTE

@test EmbeddingTE(OptimiseTraditional(), x, y, z) isa EmbeddingTE

# Internals that we may not necessarily hit using random input data.
# ----------------------------------------------------------------
emb = EmbeddingTE(OptimiseTraditional(), x, y, z)
pts, vars, τs, js = CausalityTools.te_embed(emb, x, y)
@test occursin("Tf = ", repr(vars))

# TEVars
vars1 = CausalityTools.TEVars([1], [2], [3])
vars2 = CausalityTools.TEVars([1], [2], [3], Int[])
@test vars1.C == vars2.C

# rc
# -----
x = rand(100)
X = StateSpaceSet(rand(100))

rc = CausalityTools.rc
# If multiple lags are given, and the dimension is an integer, then the number of 
# lags must match the dimension.
ds = 4
τs =  1:2
@test_throws ArgumentError rc(x, ds, τs)

# Repeated lags don't work.
ds = 2
τs =  [2, 2]
@test_throws ArgumentError rc(x, ds, τs)

# Everything should work when the dimension and the number of lags match.
ds = 2
τs =  [1, 2]
@test rc(x, ds, τs) isa Tuple # (pos, τ)

# vec of vecs input
y = [x, x]
ds = 2
τs =  [1, 2]
@test rc(y, ds, τs) isa Tuple # (pos, τ)

# different maximum dimensions per input time series
constant_τ = 3
differing_d_per_var = [2, 3]
@test rc(y, differing_d_per_var, constant_τ) isa Tuple # (pos, τ)
@test rc(y, differing_d_per_var, -constant_τ) isa Tuple # (pos, τ)

constant_τ = [-3, -2]
differing_d_per_var = [2, 3]
@test rc(y, differing_d_per_var, constant_τ, false) isa Tuple # (pos, τ)
@test rc(y, differing_d_per_var, constant_τ, true) isa Tuple # (pos, τ)

@test rc(y, ds, .-(τs)) isa Tuple # (pos, τ)

@test rc(y, 4, -1, true) isa Tuple
@test rc(y, 4, -1, false) isa Tuple