mx = 3
my = 2
mz = 4
x = rand(1:mx, 1000)
y = rand(1:my, 1000)
z = rand(2:mz+1, 1000) # ensure that starting at something other than 1 also works

c2 = contingency_matrix(x, y)
c3 = contingency_matrix(x, y, z)
@test size(c2) == (3, 2)
@test size(c3) == (3, 2, 4)
@test sum(c2) ≈ 1.0
@test sum(c3) ≈ 1.0

u = rand(1000)
v = rand(1000)
w = rand(1000)

# These are the estimators that have implementations of `marginal_encodings`
probests = [
    SymbolicPermutation(m = 3),
]

@testset "Contingency table: with $(probests[i]) discretization" for i in eachindex(probests)
    est = probests[i]
    lΩ = total_outcomes(est, u)

    # The *same* discretization is applied to all marginals. This means that the maximum
    # size of any dimension is at most `lΩ = total_outcomes(est, [x])`.
    # Not all outcomes are guaranteed to be observed, so it is possible that the size of a
    # dimension is *less* than `lΩ`.
    c2 = contingency_matrix(est, u, v)
    c3 = contingency_matrix(est, u, v, w)

    @test all(size(c2) .<= lΩ)
    @test all(size(c3) .<= lΩ)
    @test sum(c2) ≈ 1.0
    @test sum(c3) ≈ 1.0
end
