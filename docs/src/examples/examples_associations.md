# Examples of association measure estimation

## [`JointEntropyShannon`](@ref)

### [[`JointProbabilities`](@ref) with [`Dispersion`](@ref)](@id example_JointEntropyShannon_Dispersion)

```@example example_JointEntropyShannon
using CausalityTools
using Random; rng = Xoshiro(1234)
x, y = rand(rng, 100), rand(rng, 100)
measure = JointEntropyShannon()
discretization = CodifyVariables(Dispersion(m = 2, c = 3))
est = JointProbabilities(measure, discretization)
association(est, x, y)
```

## [`JointEntropyTsallis`](@ref)

### [[`JointProbabilities`](@ref) with [`OrdinalPatterns`](@ref)](@id example_JointEntropyTsallis_OrdinalPatterns)

```@example example_JointEntropyTsallis
using CausalityTools
using Random; rng = Xoshiro(1234)
x, y = rand(rng, 100), rand(rng, 100)
measure = JointEntropyTsallis()
discretization = CodifyVariables(OrdinalPatterns(m = 3))
est = JointProbabilities(measure, discretization)
association(est, x, y)
```


## [`JointEntropyRenyi`](@ref)

### [[`JointProbabilities`](@ref) with [`OrdinalPatterns`](@ref)](@id example_JointEntropyRenyi_ValueBinning)

```@example example_JointEntropyRenyi
using CausalityTools
using Random; rng = Xoshiro(1234)
x, y = rand(rng, 100), rand(rng, 100)
measure = JointEntropyRenyi(q = 0.5)
discretization = CodifyVariables(ValueBinning(2))
est = JointProbabilities(measure, discretization)
association(est, x, y)
```

## [`ConditionalEntropyShannon`](@ref)

### Analytical examples

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
disc = CodifyVariables(UniqueElements())
est = JointProbabilities(ConditionalEntropyShannon(), disc)
association(est, x, y)
```


## [`MIShannon`](@ref)

### [[`JointProbabilities`](@ref) + [`ValueBinning`](@ref)](@id example_MIShannon_JointProbabilities_ValueBinning)

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000)
discretization = CodifyVariables(ValueBinning(FixedRectangularBinning(0, 1, 5)))
est = JointProbabilities(MIShannon(), discretization)
association(est, x, y)
```


### [[`JointProbabilities`](@ref) + [`UniqueElements`](@ref)](@id example_MIShannon_JointProbabilities_UniqueElements)

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

### [Dedicated [`GaussianMI`](@ref) estimator](@id example_MIShannon_GaussianMI)

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

### [Dedicated [`KraskovStögbauerGrassberger1`](@ref) estimator](@id example_MIShannon_KSG1)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(KSG1(MIShannon(); k = 5), x, y)
```

### [Dedicated [`KraskovStögbauerGrassberger2`](@ref) estimator](@id example_MIShannon_KSG2)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(KSG2(MIShannon(); k = 5), x, y)
```

### [Dedicated [`GaoKannanOhViswanath`](@ref) estimator](@id example_MIShannon_GaoKannanOhViswanath)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(GaoKannanOhViswanath(MIShannon(); k = 10), x, y)
```

### [[`EntropyDecomposition`](@ref) + [`Kraskov`](@ref)](@id example_MIShannon_EntropyDecomposition_Kraskov)

We can compute [`MIShannon`](@ref) by naively applying a [`DifferentialEntropyEstimator`](@ref).
Note that this doesn't apply any bias correction.

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
association(EntropyDecomposition(MIShannon(), Kraskov(k = 3)), x, y)
```


### [[`EntropyDecomposition`](@ref) + [`BubbleSortSwaps`](@ref)](@id example_MIShannon_EntropyDecomposition_BubbleSortSwaps)

We can also compute [`MIShannon`](@ref) by naively applying a [`DiscreteEntropyEstimator`](@ref).
Note that this doesn't apply any bias correction.

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
disc = CodifyVariables(BubbleSortSwaps(m=5))
hest = PlugIn(Shannon())
association(EntropyDecomposition(MIShannon(), hest, disc), x, y)
```


### [[`EntropyDecomposition`](@ref) + [`Jackknife`](@ref) + [`ValueBinning`](@ref)](@id example_MIShannon_EntropyDecomposition_Jackknife_ValueBinning)

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
discretization = CodifyVariables(ValueBinning(FixedRectangularBinning(0, 1, 5)))
hest = Jackknife(Shannon())
est = EntropyDecomposition(MIShannon(), hest, discretization)
association(est, x, y)
```

## [`MIRenyiJizba`](@ref)

### [[`JointProbabilities`](@ref) + [`UniqueElements`](@ref)](@id example_MIRenyiJizba_JointProbabilities_UniqueElements)

[`MIRenyiJizba`](@ref) can be estimated for categorical data using [`JointProbabilities`](@ref) estimator
with the [`UniqueElements`](@ref) outcome space.

```@example example_mirenyijizba
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, ["a", "b", "c"], 200);
y = rand(rng, ["hello", "yoyo", "heyhey"], 200);
est = JointProbabilities(MIRenyiJizba(), UniqueElements())
association(est, x, y)
```

### [[`EntropyDecomposition`](@ref) + [`LeonenkoProzantoSavani`](@ref)](@id example_MIRenyiJizba_JointProbabilities_LeonenkoProzantoSavani)

[`MIRenyiJizba`](@ref) can also estimated for numerical data using [`EntropyDecomposition`](@ref)
in combination with any [`DifferentialInfoEstimator`](@ref) capable of estimating differential 
[`Renyi`](@ref) entropy.

```@example example_MIRenyiJizba
using CausalityTools
using Random; rng = Xoshiro(1234)
x = randn(rng, 50); y = randn(rng, 50);
def = MIRenyiJizba()
est_diff = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k=3))
association(est_diff, x, y) 
```

### [[`EntropyDecomposition`](@ref) + [`LeonenkoProzantoSavani`](@ref)](@id example_MIRenyiJizba_EntropyDecomposition_ValueBinning)

[`MIRenyiJizba`](@ref) can also estimated for numerical data using [`EntropyDecomposition`](@ref)
in combination with any [`DiscreteInfoEstimator`](@ref) capable of estimating differential 
[`Renyi`](@ref) entropy over some [`OutcomeSpace`](@ref), e.g. [`ValueBinning`](@ref).


```@example example_MIRenyiJizba
using CausalityTools
using Random; rng = Xoshiro(1234)
x = randn(rng, 50); y = randn(rng, 50);
def = MIRenyiJizba()

disc = CodifyVariables(ValueBinning(2))
est_disc = EntropyDecomposition(def, PlugIn(Renyi()), disc);
association(est_disc, x, y)
```

## [`MIRenyiSarbu`](@ref)

[`MIRenyiSarbu`](@ref) can be estimated using the [`JointProbabilities`](@ref) estimator 
in combination with any [`CodifyVariables`](@ref) or [`CodifyPoints`](@ref) discretization scheme.

### [[`JointProbabilities`](@ref) + [`UniqueElements`](@ref)](@id example_MIRenyiSarbu_JointProbabilities_UniqueElements)

```@example example_MIRenyiSarbu
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, ["a", "b", "c"], 200)
y = rand(rng, ["hello", "yoyo", "heyhey"], 200)

est = JointProbabilities(MIRenyiSarbu(), CodifyVariables(UniqueElements()))
association(est, x, y)
```

### [[`JointProbabilities`](@ref) + [`CosineSimilarityBinning`](@ref)](@id example_MIRenyiSarbu_JointProbabilities_CosineSimilarityBinning)

```@example example_MIRenyiSarbu
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 200)
y = rand(rng, 200)

est = JointProbabilities(MIRenyiSarbu(), CodifyVariables(CosineSimilarityBinning()))
association(est, x, y)
```

## [`MITsallisFuruichi`](@ref)

### [[`JointProbabilities`](@ref) + [`UniqueElements`](@ref)](@id example_MITsallisFuruichi_JointProbabilities_UniqueElements)

[`MITsallisFuruichi`](@ref) can be estimated using the [`JointProbabilities`](@ref) estimator 
in combination with any [`CodifyVariables`](@ref) or [`CodifyPoints`](@ref) discretization scheme.

```@example example_MITsallisFuruichi
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 200)
y = rand(rng, 200)

est = JointProbabilities(MITsallisFuruichi(q = 0.3), UniqueElements())
association(est, x, y) 
```

### [[`EntropyDecomposition`](@ref) + [`LeonenkoProzantoSavani`](@ref)](@id example_MITsallisFuruichi_EntropyDecomposition_LeonenkoProsantoSavani)

```@example example_MITsallisFuruichi
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 200)
y = rand(rng, 200)

est_diff = EntropyDecomposition(MITsallisFuruichi(), LeonenkoProzantoSavani(Tsallis(q= 2)))
association(est_diff, x, y)
```


### [[`EntropyDecomposition`](@ref) + [`Dispersion`](@ref)](@id example_MITsallisFuruichi_EntropyDecomposition_Dispersion)

```@example example_MITsallisFuruichi
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 200)
y = rand(rng, 200)
disc = CodifyVariables(Dispersion())
est_disc = EntropyDecomposition(MITsallisFuruichi(), PlugIn(Tsallis()), disc)

association(est_disc, x, y)
```


## [`MITsallisMartin`](@ref)

### [[`JointProbabilities`](@ref) + [`UniqueElements`](@ref)](@id example_MITsallisMartin_JointProbabilities_UniqueElements)

```@example example_MITsallisMartin
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 200)
y = rand(rng, 200)

est = JointProbabilities(MITsallisMartin(q = 1.5), UniqueElements())
association(est, x, y) 
```

### [[`EntropyDecomposition`](@ref) + [`LeonenkoProzantoSavani`](@ref)](@id example_MITsallisMartin_EntropyDecomposition_LeonenkoProsantoSavani)

[`MITsallisMartin`](@ref) can be estimated using a decomposition into entropy 
terms using [`EntropyDecomposition`](@ref) with any compatible estimator 
that can estimate differential [`Tsallis`](@ref) entropy. 


```@example example_MITsallisMartin
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 500)
y = rand(rng, 500)

est_diff = EntropyDecomposition(MITsallisMartin(), LeonenkoProzantoSavani(Tsallis(q= 1.5)))
association(est_diff, x, y)
```

### [[`EntropyDecomposition`](@ref) + [`OrdinalPatterns`](@ref)](@id example_MITsallisMartin_EntropyDecomposition_OrdinalPatterns)


```@example 
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 200)
y = rand(rng, 200)
disc = CodifyVariables(OrdinalPatterns())
est_disc = EntropyDecomposition(MITsallisMartin(), PlugIn(Tsallis()), disc)

association(est_disc, x, y)
```

## [`CMIShannon`](@ref)

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

### [[`CMIShannon`](@ref) with [`FPVP`](@ref)](@id example_CMIShannon_FPVP)

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

## [`CMIRenyiPoczos`](@ref)

### [[`PoczosSchneiderCMI`](@ref)](@id CMIRenyiPoczos_PoczosSchneiderCMI)

```@example example_cmirenyipoczos
using CausalityTools
using Distributions
using Statistics
using Random; rng = Xoshiro(1234)

n = 1000
# A chain X → Y → Z
x = rand(rng, Normal(-1, 0.5), n)
y = rand(rng, BetaPrime(0.5, 1.5), n) .+ x
z = rand(rng, Chisq(100), n)
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

## [`CMIShannon`](@ref)

### [[`MIDecomposition`](@ref) + [`KSG1`](@ref)](@id example_CMIShannon_MIDecomposition_KSG1)

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

### [[`EntropyDecomposition`](@ref) + [`Kraskov`](@ref)](@id example_CMIShannon_EntropyDecomposition_Kraskov)

Any [`DifferentialEntropyEstimator`](@ref) can also be used to compute conditional
mutual information using a sum of entropies. For that, we 
usethe [`EntropyDecomposition`](@ref) estimator. No bias correction is applied for 
[`EntropyDecomposition`](@ref) either.

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

### [[`EntropyDecomposition`](@ref) + [`ValueBinning`](@ref)](@id example_CMIShannon_EntropyDecomposition_ValueBinning)

```@example
using CausalityTools
using Distributions
n = 1000
# A chain X → Y → Z
x = rand(Epanechnikov(0.5, 1.0), n)
y = rand(Erlang(1), n) .+ x
z = rand(FDist(5, 2), n)
discretization = CodifyVariables(ValueBinning(RectangularBinning(5)))
hest = PlugIn(Shannon())
est = EntropyDecomposition(CMIShannon(), hest, discretization)
association(est, x, y, z)
```

## [`CMIRenyiJizba`](@ref)

### [[`JointProbabilities`](@ref) + [`BubbleSortSwaps`](@ref)](@id example_CMIRenyiJizba_JointProbabilities_BubbleSortSwaps)

```@example example_CMIRenyiJizba
using CausalityTools
using Random; rng = Xoshiro(1234)
x = rand(rng, 100)
y = x .+ rand(rng, 100)
z = y .+ rand(rng, 100)
disc = CodifyVariables(BubbleSortSwaps(m = 4))
est = JointProbabilities(CMIRenyiJizba(), disc)
association(est, x, z, y)
```


### [[`EntropyDecomposition`](@ref) + [`LeonenoProsantoSavani`](@ref)](@id example_CMIRenyiJizba_EntropyDecomposition_LeonenkoProzantoSavani)

```@example example_CMIRenyiJizba
using CausalityTools
using Random; rng = Xoshiro(1234)
x, y, z = rand(rng, 1000), rand(rng, 1000), rand(rng, 1000)
def = CMIRenyiJizba(q = 1.5)

# Using a differential Rényi entropy estimator
est = EntropyDecomposition(def, LeonenkoProzantoSavani(Renyi(), k = 10))
association(est, x, y, z)
```


### [[`EntropyDecomposition`](@ref) + [`OrdinalPatterns`](@ref)](@id example_CMIRenyiJizba_EntropyDecomposition_OrdinalPatterns)

```@example example_CMIRenyiJizba
using CausalityTools
using Random; rng = Xoshiro(1234)
x, y, z = rand(rng, 1000), rand(rng, 1000), rand(rng, 1000)
def = CMIRenyiJizba(q = 1.5)

# Using a plug-in Rényi entropy estimator, discretizing using ordinal patterns.
est = EntropyDecomposition(def, PlugIn(Renyi()), CodifyVariables(OrdinalPatterns(m=2)), RelativeAmount())
association(est, x, y, z)
```

## [`TEShannon`](@ref)

### [[`EntropyDecomposition`](@ref) + [`TransferOperator`](@ref)](@id example_TEShannon_EntropyDecomposition_TransferOperator)

For transfer entropy examples, we'll construct some time series for which 
there is time-delayed forcing between variables.

```@example transfer_entropy_examples

using CausalityTools
using DynamicalSystemsBase
using StableRNGs
rng = StableRNG(123)

Base.@kwdef struct Logistic4Chain{V, RX, RY, RZ, RW, C1, C2, C3, Σ1, Σ2, Σ3, RNG}
    xi::V = [0.1, 0.2, 0.3, 0.4]
    rx::RX = 3.9
    ry::RY = 3.6
    rz::RZ = 3.6
    rw::RW = 3.8
    c_xy::C1 = 0.4
    c_yz::C2 = 0.4
    c_zw::C3 = 0.35
    σ_xy::Σ1 = 0.05
    σ_yz::Σ2 = 0.05
    σ_zw::Σ3 = 0.05
    rng::RNG = Random.default_rng()
end

function eom_logistic4_chain(u, p::Logistic4Chain, t)
    (; xi, rx, ry, rz, rw, c_xy, c_yz, c_zw, σ_xy, σ_yz, σ_zw, rng) = p
    x, y, z, w = u
    f_xy = (y +  c_xy*(x + σ_xy * rand(rng)) ) / (1 + c_xy*(1+σ_xy))
    f_yz = (z +  c_yz*(y + σ_yz * rand(rng)) ) / (1 + c_yz*(1+σ_yz))
    f_zw = (w +  c_zw*(z + σ_zw * rand(rng)) ) / (1 + c_zw*(1+σ_zw))
    dx = rx * x * (1 - x)
    dy = ry * (f_xy) * (1 - f_xy)
    dz = rz * (f_yz) * (1 - f_yz)
    dw = rw * (f_zw) * (1 - f_zw)
    return SVector{4}(dx, dy, dz, dw)
end

function system(definition::Logistic4Chain)
    return DiscreteDynamicalSystem(eom_logistic4_chain, definition.xi, definition)
end

# An example system where `X → Y → Z → W`.
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 300, Ttr = 10000)))

precise = true # precise bin edges
discretization = CodifyVariables(TransferOperator(RectangularBinning(2, precise))) #
est_disc_to = EntropyDecomposition(TEShannon(), PlugIn(Shannon()), discretization);
association(est_disc_to, x, y), association(est_disc_to, y, x)
```

The Shannon-type transfer entropy from `x` to `y` is stronger than from `y` to `x`,
which is what we expect if `x` drives `y`.

```@example transfer_entropy_examples
association(est_disc_to, x, z), association(est_disc_to, x, z, y)
```

The Shannon-type transfer entropy from `x` to `z` is stronger than the transfer entropy from `x` to `z` given `y`. This is expected, because `x` drives `z` *through*
`y`, so "conditioning away" the effect of `y` should decrease the estimated 
information transfer.

## [`TERenyiJizba`](@ref)

### [[`EntropyDecomposition`](@ref) + [`TransferOperator`](@ref)](@id example_TERenyiJizba_EntropyDecomposition_TransferOperator)

We can perform the same type of analysis as above using [`TERenyiJizba`](@ref)
instead of [`TEShannon`](@ref).

```@example transfer_entropy_examples
using CausalityTools
using DynamicalSystemsBase
using StableRNGs; rng = StableRNG(123)

# An example system where `X → Y → Z → W`.
sys = system(Logistic4Chain(; rng))
x, y, z, w = columns(first(trajectory(sys, 300, Ttr = 10000)))

precise = true # precise bin edges
discretization = CodifyVariables(TransferOperator(RectangularBinning(2, precise))) #
est_disc_to = EntropyDecomposition(TERenyiJizba(), PlugIn(Renyi()), discretization);
association(est_disc_to, x, y), association(est_disc_to, y, x)
```

## [`ConvergentCrossMapping`](@ref)

We can cross map directly between two datasets.

```@example
using CausalityTools
x, y = rand(200), rand(100)
crossmap(CCM(), x, y)
```

### Cross-map with [`RandomVectors`](@ref) estimator

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

###  Cross-map with [`RandomSegments`](@ref) estimator

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
