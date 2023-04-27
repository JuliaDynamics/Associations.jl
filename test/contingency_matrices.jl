using StateSpaceSets

mx = 3
my = 2
mz = 4
x = rand(1:mx, 1000)
y = rand(1:my, 1000)
z = rand(2:mz+1, 1000) # ensure that starting at something other than 1 also works
X, Y, Z = StateSpaceSet(x), StateSpaceSet(y), StateSpaceSet(z)

# Basics
# -------------------------------------------------------------------------------------
c = contingency_matrix(x, y, z)
@test frequencies(c, dims = 1:2) isa AbstractArray{Int, 2}
@test frequencies(c, dims = 2) isa AbstractArray{Int, 1}
@test frequencies(c) isa AbstractArray{Int, 3}
@test probabilities(c, dims = 1:2) isa Probabilities{<:Real, 2}
@test probabilities(c, dims = 2) isa Probabilities{<:Real, 1}
@test probabilities(c) isa Probabilities{<:Real, 3}

c2 = contingency_matrix(x, y)
c3 = contingency_matrix(x, y, z)
@test size(c2) == (3, 2)
@test size(c3) == (3, 2, 4)
@test sum(c2) ≈ 1.0
@test sum(c3) ≈ 1.0


C = contingency_matrix(X, Y, Z)
C2 = contingency_matrix(X, Y)
@test frequencies(C, dims = 1:2) isa AbstractArray{Int, 2}
@test frequencies(C, dims = 2) isa AbstractArray{Int, 1}
@test frequencies(C) isa AbstractArray{Int, 3}
@test probabilities(C, dims = 1:2) isa Probabilities{<:Real, 2}
@test probabilities(C, dims = 2) isa Probabilities{<:Real, 1}
@test probabilities(C) isa Probabilities{<:Real, 3}
@test size(C2) == (3, 2)
@test sum(C2) ≈ 1.0

# Discretizing data before computing the contingency matrix
# -------------------------------------------------------------------------------------
u = rand(1000)
v = rand(1000)
w = rand(1000)

# These are the estimators that have implementations of `marginal_encodings`
probests = [
    SymbolicPermutation(m = 3),
    Dispersion(),
    ValueHistogram(3),
    Contingency(SymbolicPermutation(m = 3)),
    Contingency(Dispersion()),
    Contingency(ValueHistogram(3)),
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
