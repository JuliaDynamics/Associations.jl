using Test
using CausalityTools
using Random

rng = MersenneTwister(1234)
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 1000, Ttr = 10000)))

@test_throws ArgumentError asymmetry(FPVP(), x, y)
@test_throws ArgumentError estimate(PA(), x, y)
@test_throws ArgumentError estimate(PA(ηT = -1), FPVP(), x, y)
@test_throws ArgumentError estimate(PA(ηT = -1), FPVP(), x, y, z)
@test_throws ArgumentError estimate(PA(τC = 0), FPVP(), x, y)
@test_throws ArgumentError estimate(PA(τS = 0), FPVP(), x, y)
@test_throws ArgumentError estimate(PA(τC = 0), FPVP(), x, y, z)
@test_throws ArgumentError estimate(PA(τS = 0), FPVP(), x, y, z)

# Embedding optimization
p = Pecuzal(dmax = 2)
@test estimate(PA(τS = p), FPVP(), x, y) isa Vector{<:Real}
@test estimate(PA(τS = p, τC = p), FPVP(), x, y, z) isa Vector{<:Real}

# Single-value embedding parameters
@test asymmetry(PA(ηT = 1:5, τS = 1, τC = 1), FPVP(), x, y) |> length == 5
@test asymmetry(PA(ηT = 1:5, τS = 1, τC = 1), FPVP(), x, y, z) |> length == 5

α = 0.01
measure = PA(ηT = 1:5, τS = 1, τC = 1)
est = FPVP()
test = PATest(measure, est)
r_xy = independence(test, x, y)
@test r_xy isa PATestResult
@test pvalue(r_xy) < α

r_xzy = independence(test, x, z, y)
@test r_xzy isa PATestResult
@test pvalue(r_xzy) >= α


# Single-value embedding parameters
@test asymmetry(PA(ηT = 1:5, τS = [1, 2], τC = 1:2), FPVP(), x, y) |> length == 5
@test asymmetry(PA(ηT = 1:5, τS = [1, 2], τC = 1:2), FPVP(), x, y, z) |> length == 5

α = 0.01
measure = PA(ηT = 1:5, τS = [1, 2], τC = 1:2)
est = FPVP()
test = PATest(measure, est)
r_xy = independence(test, x, y)
@test r_xy isa PATestResult
@test pvalue(r_xy) < α

r_xzy = independence(test, x, z, y)
@test r_xzy isa PATestResult
@test pvalue(r_xzy) >= α
