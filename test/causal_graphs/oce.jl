using CausalityTools
using Test
using Random

rng = MersenneTwister(1234)
sys = system(Logistic4Unidir(; rng))
x, y, z, w = columns(trajectory(sys, 150, Ttr = 1000))
X = [x, y, z, w]

infer_graph(OCE(τmax = 1), X)
