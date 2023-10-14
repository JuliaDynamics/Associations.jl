# Estimating multivariate information measures


## [`MIShannon`](@ref)

### [Estimation with [`MutualInformationEstimator`](@ref)s](@id example_mi_MutualInformationEstimator)

When estimated using a [`MutualInformationEstimator`](@ref), some form of bias
correction is usually applied. The [`KraskovStögbauerGrassberger1`](@ref) and
[`KraskovStögbauerGrassberger2`](@ref) estimators are perhaps the most popular.
A common parametric estimator is [`GaussianMI`](@ref).

#### [`MIShannon`](@ref) with [`GaussianMI`](@ref)

```@example mi_demonstration
using CausalityTools
using Distributions
using Statistics

n = 1000
using CausalityTools
x = randn(1000)
y = rand(1000) .+ x
information(GaussianMI(MIShannon()), x, y) # defaults to `MIShannon()`
```

#### [`MIShannon`](@ref) with [`KraskovStögbauerGrassberger1`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
information(KSG1(MIShannon(); k = 5), x, y)
```

#### [`MIShannon`](@ref) with [`KraskovStögbauerGrassberger2`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
information(KSG2(MIShannon(); k = 5), x, y)
```

#### [`MIShannon`](@ref) with [`GaoKannanOhViswanath`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
information(GaoKannanOhViswanath(MIShannon(); k = 10), x, y)
```

#### [`MIShannon`](@ref) with [`GaoOhViswanath`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
information(GaoOhViswanath(MIShannon(); k = 10), x, y)
```

#### Reproducing Kraskov et al. (2004)

Here, we'll reproduce Figure 4 from Kraskov et al. (2004)'s seminal paper on the nearest-neighbor based mutual information estimator. We'll estimate the mutual information
between marginals of a bivariate Gaussian for a fixed time series length of 2000,
varying the number of neighbors. *Note: in the original paper, they show multiple
curves corresponding to different time series length. We only show two single curves:
one for the [`KSG1`](@ref) estimator and one for the [`KSG2`](@ref) estimator*.

```@example ex_mutualinfo
using CausalityTools
using LinearAlgebra: det
using Distributions: MvNormal
using StateSpaceSets: StateSpaceSet
using CairoMakie
using Statistics

N = 1000
c = 0.9
Σ = [1 c; c 1]
N2 = MvNormal([0, 0], Σ)
mitrue = -0.5*log(det(Σ)) # in nats
ks = [2; 5; 7; 10:10:70] .* 2

nreps = 10 # plot average over 10 independent realizations
mis_ksg1 = zeros(nreps, length(ks))
mis_ksg2 = zeros(nreps, length(ks))
for i = 1:nreps
    D2 = StateSpaceSet([rand(N2) for i = 1:N])
    X = D2[:, 1] |> StateSpaceSet
    Y = D2[:, 2] |> StateSpaceSet
    for (j, k) in enumerate(ks)
        est1 = KSG1(MIShannon(; base = ℯ); k)
        est2 = KSG2(MIShannon(; base = ℯ); k)
        mis_ksg1[i, j] = information(est1, X, Y)
        mis_ksg2[i, j] = information(est2, X, Y)
    end
end
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "k / N", ylabel = "Mutual infomation (nats)")
scatterlines!(ax, ks ./ N, mean(mis_ksg1, dims = 1) |> vec, label = "KSG1")
scatterlines!(ax, ks ./ N, mean(mis_ksg2, dims = 1) |> vec, label = "KSG2")
hlines!(ax, [mitrue], color = :black, linewidth = 3, label = "I (true)")
axislegend()
fig
```

#### [`MutualInformationEstimator`](@ref) comparison

Most estimators suffer from significant bias when applied to discrete, finite data. One possible resolution is to add a small amount of noise to discrete variables, so that the data becomes continuous in practice.

But instead of adding noise to your data, you can also consider using an
estimator that is specifically designed to deal with continuous-discrete mixture data. 
One example is the [`GaoKannanOhViswanath`](@ref) estimator. Below, we compare its
performance to [`KSG1`](@ref) on uniformly distributed discrete multivariate data.
The true mutual information is zero. While the "naive" [`KSG1`](@ref) estimator 
diverges from the true value for these data, the [`GaoKannanOhViswanath`](@ref)
converges to the true value.

```@example ex_mutualinfo
using CausalityTools
using Statistics
using StateSpaceSets: StateSpaceSet
using Statistics: mean
using CairoMakie

function compare_ksg_gkov(;
        k = 5,
        base = 2,
        nreps = 10,
        Ls = [500:100:1000; 1500; 2500; 5000; 7000])


    mis_ksg1_mix = zeros(nreps, length(Ls))
    mis_ksg1_discrete = zeros(nreps, length(Ls))
    mis_ksg1_cont = zeros(nreps, length(Ls))
    mis_gkov_mix = zeros(nreps, length(Ls))
    mis_gkov_discrete = zeros(nreps, length(Ls))
    mis_gkov_cont = zeros(nreps, length(Ls))

    for (j, L) in enumerate(Ls)
        for i = 1:nreps
            X = StateSpaceSet(float.(rand(1:8, L, 2)))
            Y = StateSpaceSet(float.(rand(1:8, L, 2)))
            Z = StateSpaceSet(rand(L, 2))
            W = StateSpaceSet(rand(L, 2))
            est_gkov = GaoKannanOhViswanath(MIShannon(; base = ℯ); k)
            est_ksg1 = KSG1(MIShannon(; base = ℯ); k)
            mis_ksg1_discrete[i, j] = mutualinfo(est_ksg1, X, Y)
            mis_gkov_discrete[i, j] = mutualinfo(est_gkov, X, Y)
            mis_ksg1_mix[i, j] = information(est_ksg1, X, Z)
            mis_gkov_mix[i, j] = information(est_gkov, X, Z)
            mis_ksg1_cont[i, j] = information(est_ksg1, Z, W)
            mis_gkov_cont[i, j] = information(est_gkov, Z, W)
        end
    end
    return mis_ksg1_mix, mis_ksg1_discrete, mis_ksg1_cont,
        mis_gkov_mix, mis_gkov_discrete, mis_gkov_cont
end

fig = Figure()
ax = Axis(fig[1, 1], 
    xlabel = "Sample size", 
    ylabel = "Mutual information (bits)")
Ls = [100; 200; 500; 1000; 2500; 5000; 7000]
nreps = 5
k = 3
mis_ksg1_mix, mis_ksg1_discrete, mis_ksg1_cont,
    mis_gkov_mix, mis_gkov_discrete, mis_gkov_cont = 
    compare_ksg_gkov(; nreps, k, Ls)

scatterlines!(ax, Ls, mean(mis_ksg1_mix, dims = 1) |> vec, 
    label = "KSG1 (mixed)", color = :black, 
    marker = :utriangle)
scatterlines!(ax, Ls, mean(mis_ksg1_discrete, dims = 1) |> vec, 
    label = "KSG1 (discrete)", color = :black, 
    linestyle = :dash, marker = '▲')
scatterlines!(ax, Ls, mean(mis_ksg1_cont, dims = 1) |> vec, 
    label = "KSG1 (continuous)", color = :black, 
    linestyle = :dot, marker = '●')
scatterlines!(ax, Ls, mean(mis_gkov_mix, dims = 1) |> vec, 
    label = "GaoKannanOhViswanath (mixed)", color = :red, 
    marker = :utriangle)
scatterlines!(ax, Ls, mean(mis_gkov_discrete, dims = 1) |> vec, 
    label = "GaoKannanOhViswanath (discrete)", color = :red, 
    linestyle = :dash, marker = '▲')
scatterlines!(ax, Ls, mean(mis_gkov_cont, dims = 1) |> vec, 
    label = "GaoKannanOhViswanath (continuous)", color = :red, 
    linestyle = :dot, marker = '●')
axislegend(position = :rb)
fig
```

### [Estimation using [`DifferentialEntropyEstimator`](@ref)s](@id example_mi_DifferentialEntropyEstimator)

#### Simple example

We can compute [`MIShannon`](@ref) by naively applying a [`DifferentialEntropyEstimator`](@ref).
Note that this doesn't apply any bias correction.

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
information(EntropyDecomposition(MIShannon(), Kraskov(k = 3)), x, y)
```

#### [`DifferentialEntropyEstimator`](@ref) comparison

Let's compare the performance of a subset of the implemented mutual information estimators. We'll use example data from Lord et al., where the analytical mutual information is known.

```@example ex_mutualinfo
using CausalityTools
using LinearAlgebra: det
using StateSpaceSets: StateSpaceSet
using Distributions: MvNormal
using LaTeXStrings
using CairoMakie

# adapted from https://juliadatascience.io/makie_colors
function new_cycle_theme()
    # https://nanx.me/ggsci/reference/pal_locuszoom.html
    my_colors = ["#D43F3AFF", "#EEA236FF", "#5CB85CFF", "#46B8DAFF",
        "#357EBDFF", "#9632B8FF", "#B8B8B8FF"]
    cycle = Cycle([:color, :linestyle, :marker], covary=true) # alltogether
    my_markers = [:circle, :rect, :utriangle, :dtriangle, :diamond,
        :pentagon, :cross, :xcross]
    my_linestyle = [nothing, :dash, :dot, :dashdot, :dashdotdot]
    return Theme(
        fontsize = 22, font="CMU Serif",
        colormap = :linear_bmy_10_95_c78_n256,
        palette = (
            color = my_colors, 
            marker = my_markers, 
            linestyle = my_linestyle,
        ),
        Axis = (
            backgroundcolor= (:white, 0.2), 
            xgridstyle = :dash, 
            ygridstyle = :dash
        ),
        Lines = (
            cycle= cycle,
        ), 
        ScatterLines = (
            cycle = cycle,
        ),
        Scatter = (
            cycle = cycle,
        ),
        Legend = (
            bgcolor = (:grey, 0.05), 
            framecolor = (:white, 0.2),
            labelsize = 13,
        )
    )
end

run(est; f::Function, # function that generates data
        base::Real = ℯ, 
        nreps::Int = 10, 
        αs = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1], 
        n::Int = 1000) =
    map(α -> information(est, f(α, n)...), αs)

function compute_results(f::Function; estimators, k = 5, k_lord = 20,
        n = 1000, base = ℯ, nreps = 10,
        as = 7:-1:0,
        αs = [1/10^(a) for a in as])
    
    is = [zeros(length(αs)) for est in estimators]
    for (k, est) in enumerate(estimators)
        tmp = zeros(length(αs))
        for i = 1:nreps
            tmp .+= run(est; f = f, αs, base, n)
        end
        is[k] .= tmp ./ nreps
    end

    return is
end

function plot_results(f::Function, ftrue::Function; 
        base, estimators, k_lord, k, 
        as = 7:-1:0, αs = [1/10^(a) for a in as], kwargs...
    )
    is = compute_results(f; 
        base, estimators, k_lord, k, as, αs, kwargs...)
    itrue = [ftrue(α; base) for α in αs]

    xmin, xmax = minimum(αs), maximum(αs)
    
    ymin = floor(Int, min(minimum(itrue), minimum(Iterators.flatten(is))))
    ymax = ceil(Int, max(maximum(itrue), maximum(Iterators.flatten(is))))
    f = Figure()
    ax = Axis(f[1, 1],
        xlabel = "α", ylabel = "I (nats)",
        xscale = log10, aspect = 1,
        xticks = (αs, [latexstring("10^{$(-a)}") for a in as]),
        yticks = (ymin:ymax)
        )
    xlims!(ax, (1/10^first(as), 1/10^last(as)))
    ylims!(ax, (ymin, ymax))
    lines!(ax, αs, itrue, 
        label = "I (true)", linewidth = 4, color = :black)
    for (i, est) in enumerate(estimators)
        if est isa EntropyDecomposition
            es = typeof(est.est).name.name |> String
        else
            es = typeof(est).name.name |> String
        end
        @show es
        lbl = occursin("Lord", es) ? "$es (k = $k_lord)" : "$es (k = $k)"
        scatter!(ax, αs, is[i], label = lbl)
        lines!(ax, αs, is[i])

    end
    axislegend()
    return f
end

set_theme!(new_cycle_theme())
k_lord = 20
k = 5
base = ℯ

def = MIShannon(base = ℯ)
estimators = [
    EntropyDecomposition(def, Kraskov(; k)),
    EntropyDecomposition(def, KozachenkoLeonenko()),
    EntropyDecomposition(def, Zhu(; k)),
    EntropyDecomposition(def, ZhuSingh(; k)),
    EntropyDecomposition(def, Gao(; k)),
    EntropyDecomposition(def, Lord(; k = k_lord)),
    EntropyDecomposition(def, LeonenkoProzantoSavani(Shannon(); k)),
    KSG1(def; k),
    KSG2(def; k),
    GaoOhViswanath(def; k),
    GaoKannanOhViswanath(def; k),
    GaussianMI(def),
]
```

#### Family 1

In this system, samples are concentrated around the diagonal $X = Y$,
and the strip of samples gets thinner as $\alpha \to 0$.

```@example ex_mutualinfo
function family1(α, n::Int)
    x = rand(n)
    v = rand(n)
    y = x + α * v
    return StateSpaceSet(x), StateSpaceSet(y)
end

# True mutual information values for these data
function ifamily1(α; base = ℯ)
    mi = -log(α) - α - log(2)
    return mi / log(base, ℯ)
end

fig = plot_results(family1, ifamily1; 
    k_lord = k_lord, k = k, nreps = 10,
    estimators = estimators,
    base = base)
```

#### Family 2

```@example ex_mutualinfo
function family2(α, n::Int)
    Σ = [1 α; α 1]
    N2 = MvNormal(zeros(2), Σ)
    D2 = StateSpaceSet([rand(N2) for i = 1:n])
    X = StateSpaceSet(D2[:, 1])
    Y = StateSpaceSet(D2[:, 2])
    return X, Y
end

function ifamily2(α; base = ℯ)
    return (-0.5 * log(1 - α^2)) / log(ℯ, base)
end

αs = 0.05:0.05:0.95
estimators = estimators
with_theme(new_cycle_theme()) do
    f = Figure();
    ax = Axis(f[1, 1], xlabel = "α", ylabel = "I (nats)")
    is_true = map(α -> ifamily2(α), αs)
    is_est = map(est -> run(est; f = family2, αs, nreps = 20), estimators)
    lines!(ax, αs, is_true, 
        label = "I (true)", color = :black, linewidth = 3)
    for (i, est) in enumerate(estimators)
        if est isa EntropyDecomposition
            estname = typeof(est.est).name.name |> String
        else
            estname = typeof(est).name.name |> String
        end
        scatterlines!(ax, αs, is_est[i], label = estname)
    end
    axislegend(position = :lt)
    return f
end
```

#### Family 3

In this system, we draw samples from a 4D Gaussian distribution distributed
as specified in the `ifamily3` function below. We let $X$ be the two first
variables, and $Y$ be the two last variables.

```@example ex_mutualinfo
function ifamily3(α; base = ℯ)
    Σ = [7 -5 -1 -3; -5 5 -1 3; -1 -1 3 -1; -3 3 -1 2+α]
    Σx = Σ[1:2, 1:2]; Σy = Σ[3:4, 3:4]
    mi = 0.5*log(det(Σx) * det(Σy) / det(Σ))
    return mi / log(ℯ, base)
end

function family3(α, n::Int)
    Σ = [7 -5 -1 -3; -5 5 -1 3; -1 -1 3 -1; -3 3 -1 2+α]
    N4 = MvNormal(zeros(4), Σ)
    D4 = StateSpaceSet([rand(N4) for i = 1:n])
    X = D4[:, 1:2]
    Y = D4[:, 3:4]
    return X, Y
end

fig = plot_results(family3, ifamily3; 
    k_lord = k_lord, k = k, nreps = 10,
    n = 2000,
    estimators = estimators, base = base)
```

We see that the [`Lord`](@ref) estimator, which estimates local volume elements using a singular-value decomposition (SVD) of local neighborhoods, outperforms the other estimators by a large margin.

### [Estimation using [`ProbabilitiesEstimator`](@ref)s](@id example_mi_ProbabilitiesEstimator)

We can also use [`ProbabilitiesEstimator`](@ref) to estimate Shannon mutual information.
This does not apply any bias correction.

#### Discrete [`MIShannon`](@ref) with [`ValueBinning`](@ref)

A [`ValueBinning`](@ref) estimator can be used to bin the data and compute
discrete Shannon mutual information.

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000)

# Use the H3-estimation method with a discrete visitation frequency based 
# probabilities estimator over a fixed grid covering the range of the data,
# which is on [0, 1].
discretization = ValueBinning(FixedRectangularBinning(0, 1, 5))
hest = PlugIn(Shannon())
est = EntropyDecomposition(MIShannon(), hest, discretization)
information(est, x, y)
```

#### Discrete [`MIShannon`](@ref) with [`JointProbabilities`](@ref) (numerical)

The above example is in fact equivalent to [`JointProbabilities`](@ref). However,
using the  [`JointProbabilities`](@ref) estimator is more flexible, because it
can also be used on [categorical data](@ref discrete_mishannon_categorical).

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000)
discretization = ValueBinning(FixedRectangularBinning(0, 1, 5))
information(JointProbabilities(MIShannon(), discretization), x, y)
```

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

discretization = UniqueElements()
est = JointProbabilities(MIShannon(), discretization)
information(est, preferences, biased_foods), information(est, preferences, random_foods)
```

## [Conditional mutual information](@id estimating_condmutualinfo)

When estimated using a [`ConditionalMutualInformationEstimator`](@ref), some form of bias
correction is usually applied. The [`FPVP`](@ref) estimator is a popular choice.

### Estimation using dedicated estimators 

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
information(FPVP(k = 5), x, z, y) # defaults to `CMIShannon()`
```

#### [`CMIShannon`](@ref) with [`MesnerShalizi`](@ref)

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
information(MesnerShalizi(; k = 10), x, z, y) # defaults to `CMIShannon()`
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
condmutualinfo(Rahimzamani(CMIShannon(base = 10); k = 10), x, z, y)
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
est = PoczosSchneiderCMI(CMIRenyiPoczos(base = 2, q = 1.2); k = 5)
condmutualinfo(est, x, z, y)
```

### Estimation using [`MIDecomposition`](@ref)

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
est = MIDecomposition(CMIShannon(base = 2), KSG1(k = 10))
condmutualinfo(est, x, z, y)
```

### Estimation using [`EntropyDecomposition`](@ref)s

Any [`DifferentialEntropyEstimator`](@ref) can also be used to compute conditional
mutual information using a sum of entropies. . For that, we 
use [`EntropyDecomposition`](@ref). No bias correction is applied for 
[`EntropyDecomposition`](@ref) either.

#### [`CMIShannon`](@ref) with [`Kraskov`](@ref)

```@example
using CausalityTools
using Distributions
n = 1000
# A chain X → Y → Z
x = rand(Epanechnikov(0.5, 1.0), n)
y = rand(Erlang(1), n) .+ x
z = rand(FDist(5, 2), n)
est = EntropyDecomposition(CMIShannon(), Kraskov(k = 5))
condmutualinfo(est, x, z, y)
```

Any [`DiscreteInfoEstimator`](@ref) that computes entropy can also be used to compute
conditional mutual information using a sum of entropies. For that, we also
use [`EntropyDecomposition`](@ref). In the discrete case, we also have to specify a
discretization (an [`OutcomeSpace`](@ref)).

#### [`CMIShannon`](@ref) with [`ValueBinning`](@ref)

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
condmutualinfo(est, x, y, z)
```
