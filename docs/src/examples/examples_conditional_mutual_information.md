# [Conditional mutual information](@id quickstart_mutualinfo)

## [`CMIShannon`](@ref)

### Estimation using [`ConditionalMutualInformationEstimator`](@ref)s

When estimated using a [`ConditionalMutualInformationEstimator`](@ref), some form of bias
correction is usually applied. The [`FPVP`](@ref) estimator is a popular choice.

#### [`CMIShannon`](@ref) with [`GaussianCMI`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
# A chain X → Y → Z
x = randn(1000)
y = randn(1000) .+ x
z = randn(1000) .+ y
condmutualinfo(GaussianCMI(), x, z, y) # defaults to `CMIShannon()`
```

#### [`CMIShannon`](@ref) with [`FPVP`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
# A chain X → Y → Z
x = rand(Normal(-1, 0.5), n)
y = rand(BetaPrime(0.5, 1.5), n) .+ x
z = rand(Chisq(100), n)
z = (z ./ std(z)) .+ y

# We expect zero (in practice: very low) CMI when computing I(X; Z | Y), because
# the link between X and Z is exclusively through Y, so when observing Y,
# X and Z should appear independent.
condmutualinfo(FPVP(k = 5), x, z, y) # defaults to `CMIShannon()`
```

#### [`CMIShannon`](@ref) with [`MesnerShalisi`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
# A chain X → Y → Z
x = rand(Normal(-1, 0.5), n)
y = rand(BetaPrime(0.5, 1.5), n) .+ x
z = rand(Chisq(100), n)
z = (z ./ std(z)) .+ y

# We expect zero (in practice: very low) CMI when computing I(X; Z | Y), because
# the link between X and Z is exclusively through Y, so when observing Y,
# X and Z should appear independent.
condmutualinfo(MesnerShalisi(k = 10), x, z, y) # defaults to `CMIShannon()`
```

#### [`CMIShannon`](@ref) with [`Rahimzamani`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
# A chain X → Y → Z
x = rand(Normal(-1, 0.5), n)
y = rand(BetaPrime(0.5, 1.5), n) .+ x
z = rand(Chisq(100), n)
z = (z ./ std(z)) .+ y

# We expect zero (in practice: very low) CMI when computing I(X; Z | Y), because
# the link between X and Z is exclusively through Y, so when observing Y,
# X and Z should appear independent.
condmutualinfo(CMIShannon(base = 10), Rahimzamani(k = 10), x, z, y)
```

#### [`CMIRenyiPoczos`](@ref) with [`PoczosSchneiderCMI`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
# A chain X → Y → Z
x = rand(Normal(-1, 0.5), n)
y = rand(BetaPrime(0.5, 1.5), n) .+ x
z = rand(Chisq(100), n)
z = (z ./ std(z)) .+ y

# We expect zero (in practice: very low) CMI when computing I(X; Z | Y), because
# the link between X and Z is exclusively through Y, so when observing Y,
# X and Z should appear independent.
condmutualinfo(CMIRenyiPoczos(base = 2, q = 1.2), PoczosSchneiderCMI(k = 5), x, z, y)
```

### Estimation using [`MutualInformationEstimator`](@ref)s

Any [`MutualInformationEstimator`](@ref) can also be used to compute conditional
mutual information using the chain rule of mutual information. However, the naive
application of these estimators don't perform any bias correction when
taking the difference of mutual information terms.

#### [`CMIShannon`](@ref) with [`KSG1`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
# A chain X → Y → Z
x = rand(Normal(-1, 0.5), n)
y = rand(BetaPrime(0.5, 1.5), n) .+ x
z = rand(Chisq(100), n)
z = (z ./ std(z)) .+ y

# We expect zero (in practice: very low) CMI when computing I(X; Z | Y), because
# the link between X and Z is exclusively through Y, so when observing Y,
# X and Z should appear independent.
condmutualinfo(CMIShannon(base = 2), KSG1(k = 5), x, z, y)
```

### Estimation using [`DifferentialEntropyEstimator`](@ref)s

Any [`DifferentialEntropyEstimator`](@ref) can also be used to compute conditional
mutual information using a sum of entropies. However, the naive
application of these estimators don't perform any bias application when
taking the sum of entropy terms.

#### [`CMIShannon`](@ref) with [`Kraskov`](@ref)

```@example
using CausalityTools
using Distributions
n = 1000
# A chain X → Y → Z
x = rand(Epanechnikov(0.5, 1.0), n)
y = rand(Erlang(1), n) .+ x
z = rand(FDist(5, 2), n)
condmutualinfo(CMIShannon(), Kraskov(k = 5), x, z, y)
```

### Estimation using [`ProbabilitiesEstimator`](@ref)s

Any [`ProbabilitiesEstimator`](@ref) can also be used to compute conditional
mutual information using a sum of entropies. However, the naive
application of these estimators don't perform any bias application when
taking the sum of entropy terms.

#### [`CMIShannon`](@ref) with [`ValueHistogram`](@ref)

```@example
using CausalityTools
using Distributions
n = 1000
# A chain X → Y → Z
x = rand(Epanechnikov(0.5, 1.0), n)
y = rand(Erlang(1), n) .+ x
z = rand(FDist(5, 2), n)
est = ValueHistogram(RectangularBinning(5))
condmutualinfo(CMIShannon(), est, x, z, y), condmutualinfo(CMIShannon(), est, x, y, z)
```
