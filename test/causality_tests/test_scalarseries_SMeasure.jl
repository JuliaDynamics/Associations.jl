
@test SMeasureTest(τ = 1) isa SMeasureTest
@test SMeasureTest(τ = 1) isa CausalityTest

# Some random time series
npts = 100
xr, yr = rand(npts), rand(npts);

# Initialise test, specifying embedding dimension, emebdding lag 
# and number of nearest neighbors
test = SMeasureTest(m = 2, τ = 1, K = 2:5)
@test causality(xr, yr, test) isa Vector{T} where T
