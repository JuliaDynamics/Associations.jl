# [Predictive asymmetry](@id examples_predictive_asymmetry)

## [Computing the asymmetry distribution](@id examples_pa_asymmetry_dist)

The following example demonstrates how to compute the [`asymmetry`](@ref) distribution
from time series input. We'll use timeseries from a chain of unidirectionally coupled
logistic maps that are coupled $X \to Y \to Z \to W$.

These examples compute the asymmetry distribution directly. Use the [`PA`](@ref)
measure with [`PATest`](@ref) for formal independence testing.

### Pairwise analysis

When considering only two variables $V_1$ and $V_2$, we expect the distribution
$\DeltaA_{X \to Y}$ to be skewed towards positive values if $V_1 \to V2$.

Parameters are tuned by providing an instance of the [`PA`](@ref)
measure, which quantifies directional influence. We'll use the [`FPVP`](@ref) estimator,
and compute the asymmetry distribution over prediction lags `ηT = 1:10`.
In real applications, it is important to ensure proper embeddings for the source
(and conditional, if relevant) variables. We will optimize embedding parameters
using the "traditional" approach from
[DelayEmbeddings.jl](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/traditional/).

```@example example_pa1
using CausalityTools
using DelayEmbeddings
using Random
rng = MersenneTwister(1234)

sys = system(Logistic4Chain(xi = [0.1, 0.2, 0.3, 0.4]; rng))
x, y, z, w = columns(trajectory(sys, 1000))
τx = estimate_delay(x, "mi_min")
τy = estimate_delay(y, "mi_min")
est = FPVP(; k = 3, w = 5)
ΔA_xy = asymmetry(PA(ηT = 1:10, τS = τx), est, x, y)
ΔA_yx = asymmetry(PA(ηT = 1:10, τS = τy), est, y, x)
ΔA_xy, ΔA_yx
```

As expected, since there is coupling $X \to Y$, $\Delta A_{X \to Y}$ is skewed
towards positive values, while $\Delta A_{Y \to X}$ is skewed towards negative values
because there is no coupling $Y \to X$.

### Conditional analysis

What happens if we compute$\Delta A_{X \to Z}$? We'd maybe expect there to be 
some information transfer $X \to Z$, even though ther are not directly linked, because
information is transferred through $Y$.

```@example example_pa1
ΔA_xz = asymmetry(PA(ηT = 1:10, τS = estimate_delay(x, "mi_min")), est, x, z)
```

As expected, the distribution is still skewed towards positive values. To determine
whether the information flow between $x$ and $z$ is mediated by $y$, we can compute
the conditional distribution $\Delta A_{X \to Z | Y}$. If these values are still positively
skewed, we conclude that $Y$ is not a mediating variable. If conditioning on $Y$ causes
$\Delta A_{X \to Z | Y}$ to not be skewed towards positive values any more, then
we conclude that $Y$ is a mediating variable and that $X$ and $Z$ are linked $X \to Y \to Z$.

In [these examples](@ref examples_patest), the same time series are formally tested
for independence using a [`PATest`](@ref).
