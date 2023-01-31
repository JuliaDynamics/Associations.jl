using Test
using Statistics

D = Dataset(rand(100, 3))
@test all(fastcov(D) .≈ cov(Matrix(D)))
@test all(fastcor(D) .≈ cor(Matrix(D)))
