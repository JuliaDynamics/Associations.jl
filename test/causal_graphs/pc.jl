using Test
using Graphs: SimpleDiGraph
using StableRNGs
rng = StableRNG(123)

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

n = 10000

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
dg_ct = infer_graph(alg, X; verbose = true)
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

dg_ct = infer_graph(alg, X; verbose = true)
dg_ci = pcalg(df, α, gausscitest)
@test dg_ct == dg_ci

# -------------------------------------------------------------------------------
# Test that different combinations of independence tests work. For this,
# we can use much shorter time series, because the purpose is just to rule
# out implementation errors, not to check that the correct result is obtained.
# -------------------------------------------------------------------------------
x, y, z = rand(rng, 50), rand(rng, 50), rand(rng, 50)
α = 0.01
X = [x, y, z]
nshuffles = 3

utests = [
    CorrTest(),
    SurrogateTest(PearsonCorrelation(); nshuffles, rng),# nonparametric version of CorrTest
    SurrogateTest(MIShannon(), KSG2(); nshuffles, rng),
    SurrogateTest(DistanceCorrelation(); nshuffles, rng),
    ];
ctests = [
    CorrTest(),
    SurrogateTest(PartialCorrelation(); nshuffles, rng), # nonparametric version of CorrTest
    LocalPermutationTest(CMIShannon(), KSG2(); nshuffles, rng),
    LocalPermutationTest(DistanceCorrelation(); nshuffles, rng),
]

tn(x) = Base.typename(typeof(x)).wrapper
for i in eachindex(combos)
    u, c = combos[i]
    @testset "PC algorithm. Pairwise: $(tn(u)). Conditional: $(tn(c))" begin
        alg = PC(u, c; α = α, maxiters_orient = 10)
        g = infer_graph(alg, X)
        @test g isa SimpleDiGraph
    end
end

alg = PC(CorrTest(), CorrTest(), maxdepth = 1)
@test infer_graph(alg, X) isa SimpleDiGraph
