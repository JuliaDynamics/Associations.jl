using Test
using Graphs: SimpleDiGraph, complete_digraph, has_edge, SimpleEdge, add_edge!
using StableRNGs
using CausalInference: pcalg, gausscitest
using Combinatorics
rng = StableRNG(123)

@test_throws ArgumentError PC(CorrTest(), CorrTest(), α = -0.5)

# -------------------------------------------------------------------------------
# "Analytical" tests
# -------------------------------------------------------------------------------
# Compare with CausalInference.jl, which already does more rigorous testing.
# Test cases 1 and 2 are taken from the documentation of CausalInference.jl.
# Tets case 3 is custom made.
# We only test using a correlation test (`CausalInference.gausscitest`),
# which we replicate here by using a `CorrTest` both for the unconditional and
# conditional case.
# -------------------------------------------------------------------------------
α = 0.01
alg = PC(CorrTest(), CorrTest(); α)

n = 1000

# Case 1
x = randn(rng, n)
v = x + randn(rng, n)*0.25
w = x + randn(rng, n)*0.25
z = v + w + randn(rng, n)*0.25
s = z + randn(rng, n)*0.25
X = [x, v, w, z, s]
df = (x=x, v=v, w=w, z=z, s=s)
dg_ct = infer_graph(alg, X; verbose = true)
dg_ci = pcalg(df, α, gausscitest)
@test dg_ct == dg_ci

# Case 2
v = randn(rng, n)
x = v + randn(rng, n)*0.25
w = x + randn(rng, n)*0.25
z = v + w + randn(rng, n)*0.25
s = z + randn(rng, n)*0.25
X = [x, v, w, z, s]
df = (x=x, v=v, w=w, z=z, s=s)
dg_ct = infer_graph(alg, X; verbose = false)
dg_ci = pcalg(df, α, gausscitest)
@test dg_ct == dg_ci

# Case 3
x = randn(rng, n)
y = x + 0.2*randn(rng, n)
z = x + 0.2*randn(rng, n)
w = y + z + 0.2*randn(rng, n)
q = w + 0.2*randn(rng, n)
r = w + 0.2*randn(rng, n)
X = [x, y, z, w, q, r]
df = (x=x, y=y, z=z, w=w, q=q, r=r)

dg_ct = infer_graph(alg, X; verbose = false)
dg_ci = pcalg(df, α, gausscitest)
@test dg_ct == dg_ci

# Case 4: a collider with marginally independent parents
n = 1000
x = randn(rng, n)              # independent root
z = randn(rng, n)              # independent root  (x ⫫ z)
y = x + z + 0.2*randn(rng, n)  # collider:  x → y ← z
X = [x, y, z]
df = (x=x, y=y, z=z)

dg_ct = infer_graph(alg, X)          # PC in this repo
dg_ci = pcalg(df, α, gausscitest)    # reference
@test dg_ct == dg_ci

# Case 5: Two independent causes a, b; i and j are common effects.
n = 1000
a = randn(rng, n)
b = randn(rng, n)
i = a .+ b .+ 0.1 .* randn(rng, n)
j = a .+ b .+ 0.1 .* randn(rng, n)
X = [a, b, i, j]
df = (a=a, b=b, i=i, j=j)

dg_ct = infer_graph(alg, X)
dg_ci = pcalg(df, α, gausscitest)
@test dg_ct == dg_ci

# ------------------------------------------
# Testing specifics steps of the algorithms 
# ------------------------------------------

# Case 1:
# --------------------------------------------------------------------------------
# `skeleton_conditional!` at level 𝓁 must test a pair whose conditioning pool
# (node i's neighbors, excluding j) has size exactly 𝓁 — the boundary case of
# Algorithm 1, line 8 in Kalisch & Bühlmann (2008), |adj(C,i)\{j}| ≥ ℓ, where
# the only available conditioning set is the full neighbor set. We call the level
# directly so the pair is evaluated at 𝓁 itself, rather than being removed at a
# lower level by the full `infer_graph` sweep.
# Data: independent a, b; common effects i = a+b, j = a+b, so i ⫫ j | {a,b}.
n = 1000
a = randn(rng, n)
b = randn(rng, n)
i = a .+ b .+ 0.1 .* randn(rng, n)
j = a .+ b .+ 0.1 .* randn(rng, n)
X = [a, b, i, j]                        # nodes: 1=a, 2=b, 3=i, 4=j

graph = complete_digraph(4)
sepset = Dict{SimpleEdge,Vector{Int}}()
Associations.skeleton_conditional!(alg, graph, sepset, X, 2)  # adj(i)\{j} = {a,b}, size 𝓁=2
@test !has_edge(graph, SimpleEdge(3, 4))   # i–j removed via i ⫫ j | {a,b}
@test !has_edge(graph, SimpleEdge(4, 3))

# Case 2: 
# --------------------------------------------------------------------------------
# `skeleton_conditional!` at level 𝓁 must test conditioning sets of size
# exactly 𝓁 (Algorithm 1, line 10 in Kalisch & Bühlmann (2008): |k| = ℓ).
# A pair separable only by a size-2 set must therefore survive 𝓁=1 and be
# removed only at 𝓁=2.
# Data: independent a, b; common effects i = a+b, j = a+b ⇒ i ⫫ j | {a,b},
# but not given a or b alone.
n = 1000
a = randn(rng, n)
b = randn(rng, n)
i = a .+ b .+ 0.1 .* randn(rng, n)
j = a .+ b .+ 0.1 .* randn(rng, n)
X = [a, b, i, j]                        # nodes: 1=a, 2=b, 3=i, 4=j

alg_wb = PC(CorrTest(), CorrTest(); α=0.01, maxdepth=Inf)   # explicit maxdepth = Inf

# 𝓁=1: adj(i)\{j} = {a,b}, so only {a} and {b} are valid conditioning sets;
# neither separates i and j, so the i–j edge must remain.
graph = complete_digraph(4)
sepset = Dict{SimpleEdge,Vector{Int}}()
Associations.skeleton_conditional!(alg_wb, graph, sepset, X, 1)
@test has_edge(graph, SimpleEdge(3, 4))
@test has_edge(graph, SimpleEdge(4, 3))

# 𝓁=2: the size-2 set {a,b} is now valid and separates i and j ⇒ edge removed.
graph = complete_digraph(4)
sepset = Dict{SimpleEdge,Vector{Int}}()
Associations.skeleton_conditional!(alg_wb, graph, sepset, X, 2)
@test !has_edge(graph, SimpleEdge(3, 4))
@test !has_edge(graph, SimpleEdge(4, 3))

# -------------------------------------------------------------------------------
# Orientation rule testing.
# -------------------------------------------------------------------------------

# Rule 3: orient X — Z into X → Z when X — Y → Z and X — W → Z
# are two chains with Y and W nonadjacent.
#   Nodes: X=1 (top), Y=2, W=3, Z=4 (bottom)
#   X—Y, X—W, X—Z undirected;  Y→Z, W→Z directed;  Y and W nonadjacent
# ---------------------------------------------------------------------------
dg = SimpleDiGraph(4)
add_edge!(dg, 1, 2);  # X — Y undirected
add_edge!(dg, 2, 1)
add_edge!(dg, 1, 3);  # X — W undirected
add_edge!(dg, 3, 1)
add_edge!(dg, 1, 4);  # X — Z undirected
add_edge!(dg, 4, 1)
add_edge!(dg, 2, 4)   # Y → Z directed
add_edge!(dg, 3, 4)   # W → Z directed
# no edge between Y=2 and W=3 ⇒ nonadjacent

changed = Associations.rule3!(alg, dg)

@test changed                                    # a rule fired
@test has_edge(dg, 1, 4)                         # X → Z retained
@test !has_edge(dg, 4, 1)                        # Z → X removed ⇒ oriented X → Z
# everything else untouched:
@test has_edge(dg, 1, 2) && has_edge(dg, 2, 1)   # X — Y still undirected
@test has_edge(dg, 1, 3) && has_edge(dg, 3, 1)   # X — W still undirected
@test has_edge(dg, 2, 4) && !has_edge(dg, 4, 2)  # Y → Z still directed
@test has_edge(dg, 3, 4) && !has_edge(dg, 4, 3)  # W → Z still directed


# -------------------------------------------------------------------------------
# Test that different combinations of independence tests work. For this,
# we can use much shorter time series, because the purpose is just to rule
# out implementation errors, not to check that the correct result is obtained.
# -------------------------------------------------------------------------------
x, y, z = rand(rng, 20), rand(rng, 20), rand(rng, 20)
α = 0.01
X = [x, y, z]
nshuffles = 2

utests = [
    CorrTest(),
    SurrogateAssociationTest(PearsonCorrelation(); nshuffles, rng),# nonparametric version of CorrTest
    SurrogateAssociationTest(KSG2(MIShannon()); nshuffles, rng),
    SurrogateAssociationTest(DistanceCorrelation(); nshuffles, rng),
    ];
ctests = [
    CorrTest(),
    SurrogateAssociationTest(PartialCorrelation(); nshuffles, rng), # nonparametric version of CorrTest
    LocalPermutationTest(MIDecomposition(CMIShannon(), KSG2()); nshuffles, rng),
    LocalPermutationTest(DistanceCorrelation(); nshuffles, rng),
]

tn(x) = Base.typename(typeof(x)).wrapper
for u in utests
    for c in ctests
        @testset "PC algorithm. Pairwise: $(tn(u)). Conditional: $(tn(c))" begin
            alg = PC(u, c; α = α, maxiters_orient = 10)
            g = infer_graph(alg, X)
            @test g isa SimpleDiGraph
        end
    end
end

alg = PC(CorrTest(), CorrTest(), maxdepth = 1)
@test infer_graph(alg, X) isa SimpleDiGraph

# In the future this should error when is_directed is implemented,
# because it shouldn't be possible to use `PC` with directed measures.
# ----------------------------------------------------------------
# x, y, z = rand(rng, 50), rand(rng, 50), rand(rng, 50)
# X = [x, y, z]
# tt = SurrogateAssociationTest(MIDecomposition(TEShannon(), KSG2()))
# ct = CorrTest()
# @test_throws ArgumentError infer_graph(PC(ct, tt), X)
# @test_throws ArgumentError infer_graph(PC(tt, ct), X)
# @test_throws ArgumentError infer_graph(PC(tt, tt), X)
