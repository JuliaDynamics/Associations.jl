# A chain of coupled logistic maps. Set the coupling from first to second 
# variable high, so that transferentropy(x → z) becomes significant.
# This should vanish when doing transferentropy(x → z | y)
sys = logistic4(c₁₂ = 0.6, u₀ = [0.1, 0.2, 0.3, 0.4]);
n = 2000
x, y,z, w = columns(trajectory(sys, n, Ttr = 10000));

α = 0.04 # Arbitrary significance level 1 - α = 0.96

# The ground truth is X → Y → Z.
test = SurrogateTest(TEShannon(), FPVP())


@test pvalue(independence(test, x, z)) < α # This has been manually tested to occur with c₁₂ = 0.8

# We should be able to reject the null when testing transferentropy(x → y | z)
@test pvalue(independence(test, x, z, y)) > α
