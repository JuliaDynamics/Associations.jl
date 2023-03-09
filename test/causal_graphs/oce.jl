using CausalityTools
using Test
using StableRNGs

rng = StableRNG(123)
sys = system(Logistic4Chain(; rng))
X = columns(trajectory(sys, 350, Ttr = 10000))

parents = infer_graph(OCE(τmax = 2), X)
@test all(x ∉ parents[1].parents_js for x in (2, 3, 4))
@test all(x ∉ parents[2].parents_js for x in (3, 4))
@test all(x ∉ parents[3].parents_js for x in (4))
@test 1 ∈ parents[2].parents_js
@test 2 ∈ parents[3].parents_js
@test 3 ∈ parents[4].parents_js
