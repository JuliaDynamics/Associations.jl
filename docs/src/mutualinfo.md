# Mutual information

## API

```@docs
mutualinfo
MutualInformationEstimator
```

## Estimators

```@docs
KraskovStögbauerGrassberger1
KraskovStögbauerGrassberger2
GaoKannanOhViswanath
GaoOhViswanath
FuruichiMI
```

## Examples

### Comparing with analytical mutual information

Let's compare the performance of a subset of the implemented mutual information estimators. We'll use example data from Lord et al., where the analytical mutual information is known.

```@example ex_mutualinfo
using CausalityTools
using LinearAlgebra: det
using StateSpaceSets: Dataset
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
    map(α -> mutualinfo(Shannon(; base), est, f(α, n)...), αs)

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
    GaoNaive(; k),
    GaoNaiveCorrected(; k),
    Lord(; k = k_lord),
    KSG1(; k), 
    KSG2(; k),
    GaoOhViswanath(; k),
    GaoKannanOhViswanath(; k),
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
    return Dataset(x), Dataset(y)
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
    D2 = Dataset([rand(N2) for i = 1:n])
    X = Dataset(D2[:, 1])
    Y = Dataset(D2[:, 2])
    return X, Y
end

function ifamily2(α; base = ℯ)
    return (-0.5 * log(1 - α^2)) / log(base, ℯ)
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
    return mi / log(base, ℯ)
end

function family3(α, n::Int)
    Σ = [7 -5 -1 -3; -5 5 -1 3; -1 -1 3 -1; -3 3 -1 2+α]
    N4 = MvNormal(zeros(4), Σ)
    D4 = Dataset([rand(N4) for i = 1:n])
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

### Reproducing Kraskov et al. (2004)

Here, we'll reproduce Figure 4 from Kraskov et al. (2004)'s seminal paper on the nearest-neighbor based mutual information estimator. We'll estimate the mutual information
between marginals of a bivariate Gaussian for a fixed time series length of 2000,
varying the number of neighbors. *Note: in the original paper, they show multiple
curves corresponding to different time series length. We only show two single curves:
one for the [`KSG1`](@ref) estimator and one for the [`KSG2`](@ref) estimator*.

```@example ex_mutualinfo
using CausalityTools
using LinearAlgebra: det
using Distributions: MvNormal
using StateSpaceSets: Dataset
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
    D2 = Dataset([rand(N2) for i = 1:N])
    X = D2[:, 1] |> Dataset
    Y = D2[:, 2] |> Dataset
    mis_ksg1[i, :] = map(k -> mutualinfo(Shannon(; base = ℯ), KSG1(; k), X, Y), ks)
    mis_ksg2[i, :] = map(k -> mutualinfo(Shannon(; base = ℯ), KSG2(; k), X, Y), ks)
end
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "k / N", ylabel = "Mutual infomation (nats)")
scatterlines!(ax, ks ./ N, mean(mis_ksg1, dims = 1) |> vec, label = "KSG1")
scatterlines!(ax, ks ./ N, mean(mis_ksg2, dims = 1) |> vec, label = "KSG2")
hlines!(ax, [mitrue], color = :black, linewidth = 3, label = "I (true)")
axislegend()
fig
```

### Continuous-discrete mixture data

Most estimators suffer from significant bias when applied to discrete
data. One possible resolution is to add a small amount of noise to discrete variables, so that the data becomes continuous in practice.

Instead of adding noise to your data, you can consider using an
estimator that is specifically designed to deal with continuous-discrete mixture data. One example is the [`GaoKannanOhViswanath`](@ref) estimator.

Here, we compare its performance to [`KSG1`](@ref) on uniformly 
distributed discrete multivariate data. The true mutual information is zero.

```@example ex_mutualinfo
using CausalityTools
using Statistics
using StateSpaceSets: Dataset
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
            X = Dataset(float.(rand(1:8, L, 2)))
            Y = Dataset(float.(rand(1:8, L, 2)))
            Z = Dataset(rand(L, 2))
            W = Dataset(rand(L, 2))
            mis_ksg1_discrete[i, j] = mutualinfo(Shannon(; base), est_ksg1, X, Y)
            mis_gkov_discrete[i, j] = mutualinfo(Shannon(; base), est_gkov, X, Y)
            mis_ksg1_mix[i, j] = mutualinfo(Shannon(; base), est_ksg1, X, Z)
            mis_gkov_mix[i, j] = mutualinfo(Shannon(; base), est_gkov, X, Z)
            mis_ksg1_cont[i, j] = mutualinfo(Shannon(; base), est_ksg1, Z, W)
            mis_gkov_cont[i, j] = mutualinfo(Shannon(; base), est_gkov, Z, W)
        end
    end
    return mis_ksg1_mix, mis_ksg1_discrete, mis_ksg1_cont,
        mis_gkov_mix, mis_gkov_discrete, mis_gkov_cont
end

fig = Figure()
ax = Axis(fig[1, 1], 
    xlabel = "Sample size", 
    ylabel = "Mutual information (bits)")
Ls = [100; 200; 500; 1000:1000:5000; 10000]
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
