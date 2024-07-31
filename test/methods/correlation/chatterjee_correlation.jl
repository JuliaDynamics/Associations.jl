using Test
using Associations
using Random
rng = Xoshiro(12345)

# ----------------------------------------------------------------
# Implementation tests
# ----------------------------------------------------------------

# We should get exactly the same results for preallocated measure 
# as for non-preallocated measure.
for i = 1:10
    x = rand(rng, 15)
    y = rand(rng, 15)

    # We must initialize identical seeds to ensure reproducible results
    rng_seed = rand(rng, 1:100) 
    m1 = ChatterjeeCorrelation(x, y, handle_ties = false, rng = Xoshiro(rng_seed))
    m2 = ChatterjeeCorrelation(handle_ties = false, rng = Xoshiro(1234))
    c1 = association(m1, x, y)
    c2 = association(m2, x, y)
    @test c1 == c2
end
# ----------------------------------------------------------------
# "analytical" test by comparing to output from the XICOR R-package by 
# the method author.
# ----------------------------------------------------------------

# Without repetitions
x = [1.2, 5.3, 2.4, 3.3, 7.7, 6.5]
y = [5.6, 6.6, 3.3, 4.4, 2.2, 15.0]
@test round(association(ChatterjeeCorrelation(x, y, handle_ties = false), x, y), digits = 6) ≈ 0.057143
@test round(association(ChatterjeeCorrelation(x, y, handle_ties = true), x, y), digits = 6) ≈ 0.057143
@test round(association(ChatterjeeCorrelation(handle_ties = false), x, y), digits = 6) ≈ 0.057143
@test round(association(ChatterjeeCorrelation(handle_ties = false), x, y), digits = 6) ≈ 0.057143

# With repetitions. We have no way of comparing directly with the XICORE
# package, because there is no option to specify the random number seed
# for splitting ties. We'll just use some data where the coefficient 
# should be positive.
x = [1.2, 5.3, 5.3, 2.4, 3.3, 7.7, 7.7, 6.5]
y = [5.6, 6.6, 3.3, 4.4, 4.4, 2.2, 2.2, 15.0] .+ x
@test round(association(ChatterjeeCorrelation(x, y, handle_ties = false), x, y), digits = 6) > 0
@test round(association(ChatterjeeCorrelation(x, y, handle_ties = true), x, y), digits = 6) > 0
@test round(association(ChatterjeeCorrelation(handle_ties = false), x, y), digits = 6) > 0
@test round(association(ChatterjeeCorrelation(handle_ties = true), x, y), digits = 6) > 0
