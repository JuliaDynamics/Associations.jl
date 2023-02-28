using CausalityTools
using Test
using Random

rng = MersenneTwister(1234)
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(trajectory(sys, 150, Ttr = 1000))
X = [x, y, z, w]

parents = infer_graph(OCE(τmax = 1), X)
@test all(x ∉ parents[1].parents_js for x in (2, 3, 4))
@test all(x ∉ parents[2].parents_js for x in (3, 4))
@test all(x ∉ parents[3].parents_js for x in (4))
@test 1 ∈ parents[2].parents_js
@test 2 ∈ parents[3].parents_js
@test 3 ∈ parents[4].parents_js
