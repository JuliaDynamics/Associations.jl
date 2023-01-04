x = rand(1:5, 1000)
y = rand(1:3, 1000)
z = rand(1:2, 1000)

c2 = ContingencyMatrix(x, y)
c3 = ContingencyMatrix(x, y, z)

@test size(c2) == (5, 3)
@test size(c3) == (5, 3, 2)
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
    c2 = ContingencyMatrix(est,u, v)
    c3 = ContingencyMatrix(est, u, v, w)

    @test all(size(c2) .<= lΩ)
    @test all(size(c3) .<= lΩ)
    @test sum(c2) ≈ 1.0
    @test sum(c3) ≈ 1.0
end
