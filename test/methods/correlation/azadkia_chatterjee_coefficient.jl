using Test
using Statistics
using Random
rng = Xoshiro(12345)
n = 15
x = rand(rng, n)
y = rand(rng, n)
z = rand(rng, n)
m = AzadkiaChatterjeeCoefficient()
@test association(m, x, y) isa Real
@test association(m, x, y, z) isa Real

# Example 8.2 from Azadkia & Chatterjee.
# We're doing 10 replicates and seeing if the mean Tₙ(Y, Z | X₁) and Tₙ(Y, Z) is within the 
# range of values they get.
Tₙs_cond = zeros(10)
Tₙs_pairwise = zeros(10)

m = AzadkiaChatterjeeCoefficient()
for i = 1:10
    n = 1000; x1, x2 = randn(n), randn(n); y = x1 .^2 .+ x2 .^ 2; z = atan.(x1 ./ x2)
    Tₙs_cond[i] = association(m, y, z, x1)
    Tₙs_pairwise[i] = association(m, y, z)
end
@test all(0.79 < mean(Tₙs_cond) < 0.84)
@test all(-0.06 < mean(Tₙs_pairwise) < 0.06)

# pretty printing
out = repr(AzadkiaChatterjeeCoefficient())
@test occursin("AzadkiaChatterjeeCoefficient", out)
@test occursin("theiler = ", out)