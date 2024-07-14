# Examples of association measure estimation

## Information measures 

### [`ConditionalEntropyShannon`](@ref): analytical example

This is essentially example 2.2.1 in Cover & Thomas (2006), where they use the following
relative frequency table as an example. Notethat Julia is column-major, so we need to
transpose their example. Then their `X` is in the first dimension of our table (along
columns) and their `Y` is our second dimension (rows).

```@example ce_contingency_table
using CausalityTools
freqs_yx = [1//8 1//16 1//32 1//32; 
    1//16 1//8  1//32 1//32;
    1//16 1//16 1//16 1//16; 
    1//4  0//1  0//1  0//1];
# `freqs_yx` is already normalized, se we can feed it directly to `Probabilities`
pxy = Probabilities(freqs_yx)
```

The marginal distribution for `x` (first dimension) is

```@example ce_contingency_table
marginal(pxy, dims = 2)
```

The marginal distribution for `y` (second dimension) is

```@example ce_contingency_table
marginal(pxy, dims = 1)
```

And the Shannon conditional entropy ``H^S(X | Y)``

```@example ce_contingency_table
ce_x_given_y = association(ConditionalEntropyShannon(), pxy) |> Rational
```

This is the same as in their example. Hooray! To compute ``H^S(Y | X)``, we just need to
flip the contingency matrix.

```@example ce_contingency_table
pyx = Probabilities(transpose(freqs_yx))
ce_y_given_x = association(ConditionalEntropyShannon(), pyx) |> Rational
```

### [`ConditionalEntropyShannon`](@ref) with [`JointProbabilities`](@ref) estimator

We can of course also estimate conditional entropy from data. To do so, we'll use the 
[`JointProbabilities`](@ref) estimator, which constructs a multivariate PMF for us.
Thus, we don't explicitly need a set of counts, like in the example above, because they
are estimated under the hood for us. 

Let's first demonstrate on some categorical data. For that, we must use
[`UniqueElements`](@ref) as the discretization (i.e. just count unique elements).

```@example
using CausalityTools, Random
rng = MersenneTwister(1234)
x = rand(rng, 1:3, 1000)
y = rand(rng, ["The Witcher", "Lord of the Rings"], 1000)
est = JointProbabilities(ConditionalEntropyShannon(), UniqueElements())
association(est, x, y)
```


### [`MIShannon`](@ref) with [`GaussianMI`](@ref) estimator

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
using CausalityTools
x = randn(1000)
y = rand(1000) .+ x
association(GaussianMI(MIShannon()), x, y) # defaults to `MIShannon()`
```

### [`MIShannon`](@ref) with [`KraskovStögbauerGrassberger1`](@ref) estimator

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(KSG1(MIShannon(); k = 5), x, y)
```

### [`MIShannon`](@ref) with [`KraskovStögbauerGrassberger2`](@ref) estimator

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(KSG2(MIShannon(); k = 5), x, y)
```

### [`MIShannon`](@ref) with [`GaoKannanOhViswanath`](@ref) estimator

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(GaoKannanOhViswanath(MIShannon(); k = 10), x, y)
```

### [`MIShannon`](@ref) with [`GaoOhViswanath`](@ref) estimator

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(GaoOhViswanath(MIShannon(); k = 10), x, y)
```

### [`MIShannon`](@ref) with [`EntropyDecomposition`](@ref) and [`DifferentialEntropyEstimator`](@ref)

We can compute [`MIShannon`](@ref) by naively applying a [`DifferentialEntropyEstimator`](@ref).
Note that this doesn't apply any bias correction.

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(EntropyDecomposition(MIShannon(), Kraskov(k = 3)), x, y)
```

### [`MIShannon`](@ref) with [`JointProbabilities`](@ref) and [`ValueBinning`](@ref)

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000)
discretization = ValueBinning(FixedRectangularBinning(0, 1, 5))
est = JointProbabilities(MIShannon(), discretization)
association(est, x, y)
```

### Discrete [`MIShannon`](@ref) with [`EntropyDecomposition`](@ref) and [`Jackknife`](@ref)

A [`ValueBinning`](@ref) estimator can be used to bin the data and compute
discrete Shannon mutual information.

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 50)
y = rand(rng, 50)

# Use the H3-estimation method with a discrete visitation frequency based 
# probabilities estimator over a fixed grid covering the range of the data,
# which is on [0, 1].
discretization = ValueBinning(FixedRectangularBinning(0, 1, 5))
hest = Jackknife(Shannon())
est = EntropyDecomposition(MIShannon(), hest, discretization)
association(est, x, y)
```

### [`MIShannon`](@ref) with [`JointProbabilities`](@ref) and [`UniqueElements`](@ref)

The [`JointProbabilities`](@ref) estimator can also be used with categorical data.
For example, let's compare the Shannon mutual information between the preferences
of a population sample with regards to different foods.

```@example mi_demonstration
using CausalityTools
n = 1000
preferences = rand(["neutral", "like it", "hate it"], n);
random_foods = rand(["water", "flour", "bananas", "booze", "potatoes", "beans", "soup"], n)
biased_foods = map(preferences) do preference
    if cmp(preference, "neutral") == 1
        return rand(["water", "flour"])
    elseif cmp(preference, "like it") == 1
        return rand(["bananas", "booze"])
    else
        return rand(["potatoes", "beans", "soup"])
    end
end

est = JointProbabilities(MIShannon(), UniqueElements())
association(est, preferences, biased_foods), association(est, preferences, random_foods)
```

### [`CMIShannon`](@ref) with [`GaussianCMI`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
# A chain X → Y → Z
x = randn(1000)
y = randn(1000) .+ x
z = randn(1000) .+ y
association(GaussianCMI(), x, z, y) # defaults to `CMIShannon()`
```

### [`CMIShannon`](@ref) with [`FPVP`](@ref)

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
association(FPVP(k = 5), x, z, y) # defaults to `CMIShannon()`
```

### [`CMIShannon`](@ref) with [`MesnerShalizi`](@ref)

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
association(MesnerShalizi(; k = 10), x, z, y) # defaults to `CMIShannon()`
```

### [`CMIShannon`](@ref) with [`Rahimzamani`](@ref)

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
association(Rahimzamani(CMIShannon(base = 10); k = 10), x, z, y)
```

### [`CMIRenyiPoczos`](@ref) with [`PoczosSchneiderCMI`](@ref)

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
est = PoczosSchneiderCMI(CMIRenyiPoczos(base = 2, q = 1.2); k = 5)
association(est, x, z, y)
```

In addition to the dedicated [`ConditionalMutualInformationEstimator`](@ref)s, any [`MutualInformationEstimator`](@ref) can also be used to compute conditional
mutual information using the chain rule of mutual information. However, the naive
application of these estimators don't perform any bias correction when
taking the difference of mutual information terms.

### [`CMIShannon`](@ref) with [`KSG1`](@ref)

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
est = MIDecomposition(CMIShannon(base = 2), KSG1(k = 10))
association(est, x, z, y)
```

Any [`DifferentialEntropyEstimator`](@ref) can also be used to compute conditional
mutual information using a sum of entropies. For that, we 
usethe [`EntropyDecomposition`](@ref) estimator. No bias correction is applied for 
[`EntropyDecomposition`](@ref) either.

### [`CMIShannon`](@ref) with [`Kraskov`](@ref)

```@example
using CausalityTools
using Distributions
n = 1000
# A chain X → Y → Z
x = rand(Epanechnikov(0.5, 1.0), n)
y = rand(Erlang(1), n) .+ x
z = rand(FDist(5, 2), n)
est = EntropyDecomposition(CMIShannon(), Kraskov(k = 5))
association(est, x, z, y)
```

Any [`DiscreteInfoEstimator`](@ref) that computes entropy can also be used to compute
conditional mutual information using a sum of entropies. For that, we also
use [`EntropyDecomposition`](@ref). In the discrete case, we also have to specify a
discretization (an [`OutcomeSpace`](@ref)).

### [`CMIShannon`](@ref) with [`ValueBinning`](@ref)

```@example
using CausalityTools
using Distributions
n = 1000
# A chain X → Y → Z
x = rand(Epanechnikov(0.5, 1.0), n)
y = rand(Erlang(1), n) .+ x
z = rand(FDist(5, 2), n)
discretization = ValueBinning(RectangularBinning(5))
hest = PlugIn(Shannon())
est = EntropyDecomposition(CMIShannon(), hest, discretization)
association(est, x, y, z)
```

## Cross-map measures

### [`ConvergentCrossMapping`](@ref) directly

```@example
using CausalityTools
x, y = rand(200), rand(100)
crossmap(CCM(), x, y)
```

### [`ConvergentCrossMapping`](@ref) with [`RandomVectors`](@ref)

When cross-mapping with the [`RandomVectors`](@ref) estimator, a single random subsample
of time indices (i.e. not in any particular order) of length `l` is drawn for each library
size `l`, and cross mapping is performed using the embedding vectors corresponding
to those time indices.

```@example
using CausalityTools
using Random; rng = MersenneTwister(1234)
x, y = randn(rng, 200), randn(rng, 200)

# We'll draw a single sample at each `l ∈ libsizes`. Sampling with replacement is then
# necessary, because our 200-pt timeseries will result in embeddings with
# less than 200 points.
est = RandomVectors(CCM(); libsizes = 50:25:200, replace = true, rng)
crossmap(est, x, y)
```

To generate a distribution of cross-map estimates for each `l ∈ libsizes`, just call
crossmap repeatedly, e.g.

```@example
using CausalityTools
using Random; rng = MersenneTwister(1234)
x, y = randn(rng, 200), randn(rng, 200)
est = RandomVectors(CCM(); libsizes = 50:25:200, replace = true, rng)
ρs = [crossmap(est, x, y) for i = 1:55]
M = hcat(ρs...)
```

Now, the `k`-th row of `M` contains `55` estimates of the correspondence measure `ρ`
at library size `libsizes[k]`.

### [`ConvergentCrossMapping`](@ref) with [`RandomSegments`](@ref)

When cross-mapping with the [`RandomSegments`](@ref) estimator, a single random subsample
of continguous, ordered time indices of length `l` is drawn for each library
size `l`, and cross mapping is performed using the embedding vectors corresponding
to those time indices.

```@example
using CausalityTools
using Random; rng = MersenneTwister(1234)
x, y = randn(rng, 200), randn(rng, 200)

# We'll draw a single sample at each `l ∈ libsizes`. We limit the library size to 100, 
# because drawing segments of the data longer than half the available data doesn't make
# much sense.
est = RandomSegment(CCM(); libsizes = 50:25:100, rng)
crossmap(est, x, y)
```

As above, to generate a distribution of cross-map estimates for each `l ∈ libsizes`, just call
crossmap repeatedly, e.g.

```@example
using CausalityTools
using Random; rng = MersenneTwister(1234)
x, y = randn(rng, 200), randn(rng, 200)
est = RandomSegment(CCM(); libsizes = 50:25:100, rng)
ρs = [crossmap(est, x, y) for i = 1:80]
M = hcat(ρs...)
```

Now, the `k`-th row of `M` contains `80` estimates of the correspondence measure `ρ`
at library size `libsizes[k]`.
