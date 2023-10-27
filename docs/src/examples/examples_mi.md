# [Mutual information](@id quickstart_mutualinfo)

## [`MIShannon`](@ref)

### [Estimation using [`MutualInformationEstimator`](@ref)s](@id example_mi_MutualInformationEstimator)

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
mutualinfo(KSG1(k = 5), x, y)
mutualinfo(GaussianMI(), x, y) # defaults to `MIShannon()`
```

#### [`MIShannon`](@ref) with [`KraskovStögbauerGrassberger1`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
mutualinfo(KSG1(k = 5), x, y)
```

#### [`MIShannon`](@ref) with [`KraskovStögbauerGrassberger2`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
mutualinfo(KSG2(k = 5), x, y)
```

#### [`MIShannon`](@ref) with [`GaoKannanOhViswanath`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
mutualinfo(GaoKannanOhViswanath(k = 10), x, y)
```

#### [`MIShannon`](@ref) with [`GaoOhViswanath`](@ref)

```@example mi_demonstration
using CausalityTools
x, y = rand(1000), rand(1000)
mutualinfo(GaoOhViswanath(k = 10), x, y)
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

N = 2000
c = 0.9
Σ = [1 c; c 1]
N2 = MvNormal([0, 0], Σ)
mitrue = -0.5*log(det(Σ)) # in nats
ks = [2; 5; 7; 10:10:70] .* 2

nreps = 30
mis_ksg1 = zeros(nreps, length(ks))
mis_ksg2 = zeros(nreps, length(ks))
for i = 1:nreps
    D2 = StateSpaceSet([rand(N2) for i = 1:N])
    X = D2[:, 1] |> StateSpaceSet
    Y = D2[:, 2] |> StateSpaceSet
    measure = MIShannon(; base = ℯ)
    mis_ksg1[i, :] = map(k -> mutualinfo(measure, KSG1(; k), X, Y), ks)
    mis_ksg2[i, :] = map(k -> mutualinfo(measure, KSG2(; k), X, Y), ks)
end
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "k / N", ylabel = "Mutual information (nats)")
scatterlines!(ax, ks ./ N, mean(mis_ksg1, dims = 1) |> vec, label = "KSG1")
scatterlines!(ax, ks ./ N, mean(mis_ksg2, dims = 1) |> vec, label = "KSG2")
hlines!(ax, [mitrue], color = :black, linewidth = 3, label = "I (true)")
axislegend()
fig
```

#### [`MutualInformationEstimator`](@ref) comparison

Most estimators suffer from significant bias when applied to discrete, finite data. One possible resolution is to add a small amount of noise to discrete variables, so that the data becomes continuous in practice.

Instead of adding noise to your data, you can consider using an
estimator that is specifically designed to deal with continuous-discrete mixture data. One example is the [`GaoKannanOhViswanath`](@ref) estimator.

Here, we compare its performance to [`KSG1`](@ref) on uniformly
distributed discrete multivariate data. The true mutual information is zero.

```@example ex_mutualinfo
using CausalityTools
using Statistics
using StateSpaceSets: StateSpaceSet
using Statistics: mean
using CairoMakie

function compare_ksg_gkov(;
        k = 5,
        base = 2,
        nreps = 15,
        Ls = [500:100:1000; 1500; 2000; 3000; 4000; 5000; 1000])

    est_gkov = GaoKannanOhViswanath(; k)
    est_ksg1 = KSG1(; k)

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
            measure = MIShannon(; base = ℯ)
            mis_ksg1_discrete[i, j] = mutualinfo(measure, est_ksg1, X, Y)
            mis_gkov_discrete[i, j] = mutualinfo(measure, est_gkov, X, Y)
            mis_ksg1_mix[i, j] = mutualinfo(measure, est_ksg1, X, Z)
            mis_gkov_mix[i, j] = mutualinfo(measure, est_gkov, X, Z)
            mis_ksg1_cont[i, j] = mutualinfo(measure, est_ksg1, Z, W)
            mis_gkov_cont[i, j] = mutualinfo(measure, est_gkov, Z, W)
        end
    end
    return mis_ksg1_mix, mis_ksg1_discrete, mis_ksg1_cont,
        mis_gkov_mix, mis_gkov_discrete, mis_gkov_cont
end

fig = Figure()
ax = Axis(fig[1, 1], 
    xlabel = "Sample size", 
    ylabel = "Mutual information (bits)")
Ls = [100; 200; 500; 1000; 2500; 5000; 10000]
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
mutualinfo(Kraskov(k = 3), x, y)
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
    cycle = Cycle([:color, :linestyle, :marker], covary=true) # all together
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
    map(α -> mutualinfo(MIShannon(; base), est, f(α, n)...), αs)

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
        es = string(typeof(est).name.name)
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

estimators = [
    Kraskov(; k), 
    KozachenkoLeonenko(),
    Zhu(; k), 
    ZhuSingh(; k),
    Gao(; k),
    Lord(; k = k_lord),
    KSG1(; k), 
    KSG2(; k),
    GaoOhViswanath(; k),
    GaoKannanOhViswanath(; k),
    GaussianMI(),
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
        estname = typeof(est).name.name |> String
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

#### Discrete [`MIShannon`](@ref) with [`ValueHistogram`](@ref)

A [`ValueHistogram`](@ref) estimator can be used to bin the data and compute
discrete Shannon mutual information.

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000)

# Use the H3-estimation method with a discrete visitation frequency based 
# probabilities estimator over a fixed grid covering the range of the data,
# which is on [0, 1].
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(est, x, y)
```

#### Discrete [`MIShannon`](@ref) with [`Contingency`](@ref) (numerical)

The above example is in fact equivalent to [`Contingency`](@ref). However,
using the  [`Contingency`](@ref) estimator is more flexible, because it
can also be used on [categorical data](@ref discrete_mishannon_categorical).

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000)
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(Contingency(est), x, y)
```

#### Discrete [`MIShannon`](@ref) with [`ContingencyMatrix`](@ref) (manual)

If you need explicit access to the estimated joint probability mass function,
use a [`ContingencyMatrix`](@ref) directly.

```@example mi_demonstration
using CausalityTools
using Random; rng = MersenneTwister(1234)
x = rand(rng, 1000)
y = rand(rng, 1000)
c = contingency_matrix(est, x, y)
est = ValueHistogram(FixedRectangularBinning(0, 1, 5))
mutualinfo(c)
```

#### [Discrete [`MIShannon`](@ref) with [`Contingency`](@ref) (categorical)](@id discrete_mishannon_categorical)

The [`ContingencyMatrix`](@ref) approach can also be used with categorical data.
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

c_biased = contingency_matrix(preferences, biased_foods) 
c_random = contingency_matrix(preferences, random_foods) 
mutualinfo(c_biased), mutualinfo(c_random)
```

#### Longer example: AR1-system and unidirectionally coupled logistic maps

In this example we generate realizations of two different systems where we know the strength of coupling between the variables. Our aim is to compute Shannon mutual information $I^S(X; Y)$ ([`MIShannon`](@ref)) between time series of each variable and assess how the magnitude of $I^S(X; Y)$ changes as we change the strength of coupling between $X$ and $Y$. We'll use two systems that ship with CausalityTools.jl:

* A stochastic system consisting of two unidirectionally coupled first-order autoregressive processes ([`ar1_unidir`](@ref))
* A deterministic, chaotic system consisting of two unidirectionally coupled logistic maps ([`logistic2_unidir`](@ref))

We use the default input parameter values (see [`AR1Unidir`](@ref) and [`Logistic2Unidir`](@ref) for details) and below we toggle only the random initial conditions and the coupling strength parameter `c_xy`. For each value of `c_xy` we generate 1,000 unique realizations of the system and obtain 500-point time series of the coupled variables.

To estimate the mutual information, we use the binning-based [`ValueHistogram`](@ref) estimator. We summarize the distribution of $I(X; Y)$ values across all realizations using the median and quantiles encompassing 95 % of the values.

```@example
using CausalityTools
using Statistics
using CairoMakie

# Span a range of x-y coupling strengths
c = 0.0:0.1:1.0

# Number of observations in each time series
npts = 500

# Number of unique realizations of each system
n_realizations = 1000

# Get MI for multiple realizations of two systems, 
# saving three quantiles for each c value
mi = zeros(length(c), 3, 2)

# Define an estimator for MI
b = RectangularBinning(4)
estimator = ValueHistogram(b)

for i in 1 : length(c)
    
    tmp = zeros(n_realizations, 2)
    
    for k in 1 : n_realizations
        
        # Obtain time series realizations of the two 2D systems 
        # for a given coupling strength and random initial conditions
        s_logistic = system(Logistic2Unidir(; xi = rand(2), c_xy = c[i]))
        s_ar = system(AR1Unidir(xi = rand(2), c_xy = c[i]))
        lmap = first(trajectory(s_logistic, npts - 1, Ttr = 500))
        ar1 = first(trajectory(s_ar, npts - 1))
        
        # Compute the MI between the two coupled components of each system
        tmp[k, 1] = mutualinfo(MIShannon(), estimator, lmap[:, 1], lmap[:, 2])
        tmp[k, 2] = mutualinfo(MIShannon(), estimator, ar1[:, 1], ar1[:, 2])
    end
    
    # Compute lower, middle, and upper quantiles of MI for each coupling strength
    mi[i, :, 1] = quantile(tmp[:, 1], [0.025, 0.5, 0.975])
    mi[i, :, 2] = quantile(tmp[:, 2], [0.025, 0.5, 0.975])
end

# Plot distribution of MI values as a function of coupling strength for both systems
fig = with_theme(theme_minimal()) do
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "Coupling strength", ylabel = "Mutual information")
    band!(ax, c, mi[:, 1, 1], mi[:, 3, 1], color = (:black, 0.3))
    lines!(ax, c, mi[:, 2, 1], label = "2D chaotic logistic maps", color = :black)
    band!(ax, c, mi[:, 1, 2], mi[:, 3, 2], color = (:red, 0.3))
    lines!(ax, c, mi[:, 2, 2],  label = "2D order-1 autoregressive", color = :red)
    return fig
end
fig
```

As expected, $I(X; Y)$ increases with coupling strength in a system-specific manner.
