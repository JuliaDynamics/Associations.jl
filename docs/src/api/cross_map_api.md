```@meta
CollapsedDocStrings = true
```

# Cross-map API

## Estimators

```@docs
CrossmapEstimator
RandomVectors
RandomSegment
ExpandingSegment
Ensemble
```

## Advanced utility methods

For most use cases, it is sufficient to provide a [`CrossmapEstimator`](@ref) to 
[`association`](@ref) to compute a cross map measure. However, in some cases it 
can be useful to have more fine-grained controls. We offer a few utility functions
for this purpose.

These functions are used in the examples below, where we [reproduce Figures 3C and 3D](@ref example_sugihara_figs3Cand3D) of [Sugihara2012](@citet) and [reproduce figures](@ref example_PairwiseAsymmetricInference_reproduce_mccracken) from [McCracken2014](@citet).
```@docs
predict
crossmap
```


### [Example: reproducing Sugihara et al. (2012)](@id example_ConvergentCrossMapping_reproducing_sugihara)

!!! note "Run blocks consecutively"
    If copying these examples and running them locally, make sure the relevant packages (given in the first block) are loaded first.

#### Figure 3A

Let's reproduce figure 3A too, focusing only on [`ConvergentCrossMapping`](@ref) this time. In this figure, they compute the cross mapping for libraries of increasing size, always starting at time index 1. This approach - which we here call the [`ExpandingSegment`](@ref) estimator - is one of many ways of estimating the correspondence between observed and predicted value.

For this example, they use a bidirectional system with asymmetrical coupling strength.

```@example MAIN_CCM
using CausalityTools
using Statistics
using LabelledArrays
using StaticArrays
using DynamicalSystemsBase
using StateSpaceSets
using CairoMakie, Printf

function eom_logistic_sugi(u, p, t)
    (; rx, ry, βxy, βyx) = p
    (; x, y) = u

    dx = x*(rx - rx*x - βxy*y)
    dy = y*(ry - ry*y - βyx*x)
    return SVector{2}(dx, dy)
end

# βxy := effect on x of y
# βyx := effect on y of x
function logistic_sugi(; u0 = rand(2), rx, ry, βxy, βyx)
    p = @LArray [rx, ry, βxy, βyx] (:rx, :ry, :βxy, :βyx)
    DiscreteDynamicalSystem(eom_logistic_sugi, u0, p)
end

# Used in `reproduce_figure_3A_naive`, and `reproduce_figure_3A_ensemble` below.
function add_to_fig!(fig_pos, libsizes, ρs_x̂y, ρs_ŷx; title = "", quantiles = false)
    ax = Axis(fig_pos; title, aspect = 1,
        xlabel = "Library size", ylabel = "Correlation (ρ)")
    ylims!(ax, (-1, 1))
    hlines!([0], linestyle = :dash, alpha = 0.5, color = :grey)
    scatterlines!(libsizes, median.(ρs_x̂y), label = "x̂|y", color = :blue)
    scatterlines!(libsizes, median.(ρs_ŷx), label = "ŷ|x", color = :red)
    if quantiles
        band!(libsizes, quantile.(ρs_x̂y, 0.05), quantile.(ρs_x̂y, 0.95), color = (:blue, 0.5))
        band!(libsizes, quantile.(ρs_ŷx, 0.05), quantile.(ρs_ŷx, 0.95), color = (:red, 0.5))
    end
    axislegend(ax, position = :rb)
end

function reproduce_figure_3A_naive(definition::CrossmapMeasure)
    sys_bidir = logistic_sugi(; u0 = [0.2, 0.4], rx = 3.7, ry = 3.700001, βxy = 0.02, βyx = 0.32);
    x, y = columns(first(trajectory(sys_bidir, 3100, Ttr = 10000)));
    libsizes = [20:2:50; 55:5:200; 300:50:500; 600:100:900; 1000:500:3000]
    est = ExpandingSegment(definition; libsizes);
    ρs_x̂y = crossmap(est, x, y)
    ρs_ŷx = crossmap(est, y, x)

    with_theme(theme_minimal(),
        markersize = 5) do
        fig = Figure(resolution = (800, 300))
        add_to_fig!(fig[1, 1], libsizes, ρs_x̂y, ρs_ŷx; title = "`ExpandingSegment`")
        fig
    end
end

reproduce_figure_3A_naive(ConvergentCrossMapping(d = 3))
```

Hm. This looks a bit like the paper, but the curve is not smooth. We can do better!

It is not clear from the paper exactly *what* they plot in their Figure 3A, if they plot an average of some kind, or precisely what parameters and initial conditions they use. However, we can get a smoother plot by using a [`Ensemble`](@ref). Combined with a [`CrossmapEstimator`](@ref), it uses Monte Carlo resampling on subsets of the input data to compute an ensemble of `ρ`s that we here use to compute the median and 90-th percentile range for each library size.

```@example MAIN_CCM
function reproduce_figure_3A_ensemble(definition::CrossmapMeasure)
    sys_bidir = logistic_sugi(; u0 = [0.4, 0.2], rx = 3.8, ry = 3.5, βxy = 0.02, βyx = 0.1);
    x, y = columns(first(trajectory(sys_bidir, 5000, Ttr = 10000)));
    # Note: our time series are 1000 points long. When embedding, some points are
    # lost, so we must use slightly less points for the segments than 
    # there are points in the original time series.
    libsizes = [20:5:50; 55:5:200; 300:50:500; 600:100:900; 1000:500:2000]
    # No point in doing more than one rep, because there data are always the same
    # for `ExpandingSegment.`
    ensemble_ev = Ensemble(ExpandingSegment(definition; libsizes); nreps = 1)
    ensemble_rs = Ensemble(RandomSegment(definition; libsizes); nreps = 30)
    ensemble_rv = Ensemble(RandomVectors(definition; libsizes); nreps = 30)
    ρs_x̂y_es = crossmap(ensemble_ev, x, y)
    ρs_ŷx_es = crossmap(ensemble_ev, y, x)
    ρs_x̂y_rs = crossmap(ensemble_rs, x, y)
    ρs_ŷx_rs = crossmap(ensemble_rs, y, x)
    ρs_x̂y_rv = crossmap(ensemble_rv, x, y)
    ρs_ŷx_rv = crossmap(ensemble_rv, y, x)

    with_theme(theme_minimal(),
        markersize = 5) do
        fig = Figure(resolution = (800, 300))
        add_to_fig!(fig[1, 1], libsizes, ρs_x̂y_es, ρs_ŷx_es; title = "`ExpandingSegment`", quantiles = false) # quantiles make no sense for `ExpandingSegment`
        add_to_fig!(fig[1, 2], libsizes, ρs_x̂y_rs, ρs_ŷx_rs; title = "`RandomSegment`", quantiles = true)
        add_to_fig!(fig[1, 3], libsizes, ρs_x̂y_rv, ρs_ŷx_rv; title = "`RandomVector`", quantiles = true)
        fig
    end
end

reproduce_figure_3A_ensemble(ConvergentCrossMapping(d = 3, τ = -1))
```

With the [`RandomVectors`](@ref) estimator, the mean of our ensemble `ρ`s seem to look pretty much identical to Figure 3A in Sugihara et al. The [`RandomSegment`](@ref) estimator also performs pretty well, but since subsampled segments are contiguous, there are probably some autocorrelation effects at play.

We can avoid the autocorrelation issue by tuning the `w` parameter of the [`ConvergentCrossMapping`](@ref) measure, which is the 
[Theiler window](https://juliadynamics.github.io/DynamicalSystems.jl/dev/embedding/StateSpaceSet/#Theiler-window). Setting the Theiler window to `w > 0`, we can exclude neighbors of a query point `p` that are close to `p` in time, and thus deal with autocorrelation issues that way (the default `w = 0` excludes only the point itself). Let's re-do the analysis with `w = 5`, just for fun.

```@example MAIN_CCM
reproduce_figure_3A_ensemble(ConvergentCrossMapping(d = 3, τ = -1, w = 5))
```

There wasn't really that much of a difference, since for the logistic map, the autocorrelation function flips sign for every lag increase. However, for examples from other systems, tuning `w` may be important.


#### Figure 3B

What about figure 3B? Here they generate time series of length 400 for a range of values for both coupling parameters, and plot the dominant direction $\Delta = \rho(\hat{x} | y) - \rho(\hat{y} | x)$.

In the paper, they use a 1000 different parameterizations for the logistic map parameters, but don't state what is summarized in the plot. For simplicity, we'll therefore just stick to `rx = ry = 3.7`, as in the examples above, and just loop over the coupling strengths in either direction.

```@example MAIN_CCM
function reproduce_figure_3B()
    βxys = 0.0:0.02:0.4
    βyxs = 0.0:0.02:0.4
    ρx̂ys = zeros(length(βxys), length(βyxs))
    ρŷxs = zeros(length(βxys), length(βyxs))

    for (i, βxy) in enumerate(βxys)
        for (j, βyx) in enumerate(βyxs)
            sys_bidir = logistic_sugi(; u0 = [0.2, 0.4], rx = 3.7, ry = 3.7, βxy, βyx);
            # Generate 1000 points. Randomly select a 400-pt long segment.
            x, y = columns(first(trajectory(sys_bidir, 400, Ttr = 10000)));
            definition = CCM(d = 3, w = 5, τ = -1)
            ensemble = Ensemble(RandomVectors(definition; libsizes = 100), nreps = 50)
            ρx̂ys[i, j] = mean(crossmap(ensemble, x, y))
            ρŷxs[i, j] = mean(crossmap(ensemble, y, x))
        end
    end
    Δ = ρŷxs .- ρx̂ys

    with_theme(theme_minimal(),
        markersize = 5) do
        fig = Figure();
        ax = Axis(fig[1, 1], xlabel = "βxy", ylabel = "βyx")
        cont = contourf!(ax, Δ, levels = range(-1, 1, length = 10),
            colormap = :curl)
        ax.xticks = 1:length(βxys), string.([i % 2 == 0 ? βxys[i] : "" for i in 1:length(βxys)])
        ax.yticks = 1:length(βyxs), string.([i % 2 == 0 ? βyxs[i] : "" for i in 1:length(βyxs)])
        Colorbar(fig[1 ,2], cont, label = "Δ (ρ(ŷ|x) - ρ(x̂|y))")
        tightlimits!(ax)
        fig
    end
end

reproduce_figure_3B()
```

#### [Figures 3C and 3D](@id example_sugihara_figs3Cand3D)

Let's reproduce figures 3C and 3D in [Sugihara2012](@citet), which
introduced the [`ConvergentCrossMapping`](@ref) measure.
Equations and parameters can be found in their supplementary material.
Simulatenously, we also compute the [`PairwiseAsymmetricInference`](@ref) measure
from [McCracken2014](@citet), which is a related method, but uses a
slightly different embedding.


```@example MAIN_CCM
using CausalityTools
using Statistics
using LabelledArrays
using StaticArrays
using DynamicalSystemsBase
using StateSpaceSets
using CairoMakie, Printf

# -----------------------------------------------------------------------------------------
# Create 500-point long time series for Sugihara et al. (2012)'s example for figure 3.
# -----------------------------------------------------------------------------------------
sys_unidir = logistic_sugi(; u0 = [0.2, 0.4], rx = 3.7, ry = 3.700001, βxy = 0.00, βyx = 0.32);
x, y = columns(first(trajectory(sys_unidir, 500, Ttr = 10000)));

# -----------------------------------------------------------------------------------------
# Cross map.
# -----------------------------------------------------------------------------------------
m_ccm = ConvergentCrossMapping(d = 2)
m_pai = PairwiseAsymmetricInference(d = 2)
# Make predictions x̂y, i.e. predictions `x̂` made from embedding of y (AND x, if PAI)
t̂ccm_x̂y, tccm_x̂y, ρccm_x̂y = predict(m_ccm, x, y)
t̂pai_x̂y, tpai_x̂y, ρpai_x̂y = predict(m_pai, x, y);
# Make predictions ŷx, i.e. predictions `ŷ` made from embedding of x (AND y, if PAI)
t̂ccm_ŷx, tccm_ŷx, ρccm_ŷx = predict(m_ccm, y, x)
t̂pai_ŷx, tpai_ŷx, ρpai_ŷx = predict(m_pai, y, x);

# -----------------------------------------------------------------------------------------
# Plot results
# -----------------------------------------------------------------------------------------
ρs = (ρccm_x̂y, ρpai_x̂y, ρccm_ŷx, ρpai_ŷx)
sccm_x̂y, spai_x̂y, sccm_ŷx, spai_ŷx = (map(ρ -> (@sprintf "%.3f" ρ), ρs)...,)

ρs = (ρccm_x̂y, ρpai_x̂y, ρccm_ŷx, ρpai_ŷx)
sccm_x̂y, spai_x̂y, sccm_ŷx, spai_ŷx = (map(ρ -> (@sprintf "%.3f" ρ), ρs)...,)

with_theme(theme_minimal(),
    markersize = 5) do
    fig = Figure();
    ax_ŷx = Axis(fig[2,1], aspect = 1, xlabel = "y(t) (observed)", ylabel = "ŷ(t) | x (predicted)")
    ax_x̂y = Axis(fig[2,2], aspect = 1, xlabel = "x(t) (observed)", ylabel = "x̂(t) | y (predicted)")
    xlims!(ax_ŷx, (0, 1)), ylims!(ax_ŷx, (0, 1))
    xlims!(ax_x̂y, (0, 1)), ylims!(ax_x̂y, (0, 1))
    ax_ts = Axis(fig[1, 1:2], xlabel = "Time (t)", ylabel = "Value")
    scatterlines!(ax_ts, x[1:300], label = "x")
    scatterlines!(ax_ts, y[1:300], label = "y")
    axislegend()
    scatter!(ax_ŷx, tccm_ŷx, t̂ccm_ŷx, label = "CCM (ρ = $sccm_ŷx)", color = :black)
    scatter!(ax_ŷx, tpai_ŷx, t̂pai_ŷx, label = "PAI (ρ = $spai_ŷx)", color = :red)
    axislegend(ax_ŷx, position = :lt)
    scatter!(ax_x̂y, tccm_x̂y, t̂ccm_x̂y, label = "CCM (ρ = $sccm_x̂y)", color = :black)
    scatter!(ax_x̂y, tpai_x̂y, t̂pai_x̂y, label = "PAI (ρ = $spai_x̂y)", color = :red)
    axislegend(ax_x̂y, position = :lt)
    fig
end
```


### [Example: reproducing McCracken & Weigel (2014)](@id example_PairwiseAsymmetricInference_reproduce_mccracken)

Let's try to reproduce figure 8 from [McCracken2014](@citet)'s
paper on [`PairwiseAsymmetricInference`](@ref) (PAI). We'll start by defining the their example B (equations 6-7). This system consists of two
variables ``X`` and ``Y``, where ``X`` drives ``Y``.

After we have computed the PAI in both directions, we define a measure of directionality as the difference between PAI in the ``X \to Y`` direction and in the ``Y \to X`` direction, so that if ``X`` drives ``Y``, then ``\Delta < 0``.

```@example MAIN_CCM
using CausalityTools
using LabelledArrays
using StaticArrays
using DynamicalSystemsBase
using StateSpaceSets
using CairoMakie, Printf
using Distributions: Normal
using Statistics: mean, std

function eom_nonlinear_sindriver(dx, x, p, n)
    a, b, c, t, Δt = (p...,)
    x, y = x[1], x[2]
    𝒩 = Normal(0, 1)
    
    dx[1] = sin(t)
    dx[2] = a*x * (1 - b*x) + c* rand(𝒩)
    p[end-1] += 1 # update t

    return
end

function nonlinear_sindriver(;u₀ = rand(2), a = 1.0, b = 1.0, c = 2.0, Δt = 1)
    DiscreteDynamicalSystem(eom_nonlinear_sindriver, u₀, [a, b, c, 0, Δt])
end

function reproduce_figure_8_mccraken(; 
        c = 2.0, Δt = 0.2,
        as = 0.5:0.5:5.0,
        bs = 0.5:0.5:5.0)
    # -----------------------------------------------------------------------------------------
    # Generate many time series for many different values of the parameters `a` and `b`,
    # and compute PAI. This will replicate the upper right panel of 
    # figure 8 in McCracken & Weigel (2014).
    # -----------------------------------------------------------------------------------------
    
    measure = PairwiseAsymmetricInference(d = 3)

    # Manually resample `nreps` length-`L` time series and use mean ρ(x̂|X̄y) - ρ(ŷ|Ȳx)
    # for each parameter combination.
    nreps = 50
    L = 200 # length of timeseries
    Δ = zeros(length(as), length(bs))
    for (i, a) in enumerate(as)
        for (j, b) in enumerate(bs)
            s = nonlinear_sindriver(; a, b, c,  Δt)
            x, y = columns(first(trajectory(s, 1000, Ttr = 10000)))
            Δreps = zeros(nreps)
            for i = 1:nreps
                # Ensure we're subsampling at the same time indices. 
                ind_start = rand(1:(1000-L))
                r = ind_start:(ind_start + L)
                Δreps[i] = @views crossmap(measure, y[r], x[r]) - 
                    crossmap(measure, x[r], y[r])
            end
            Δ[i, j] = mean(Δreps)
        end
    end

    # -----------------------------------------------------------------------------------------
    # An example time series for plotting.
    # -----------------------------------------------------------------------------------------
    sys = nonlinear_sindriver(; a = 1.0, b = 1.0, c, Δt)
    npts = 500
    orbit = first(trajectory(sys, npts, Ttr = 10000))
    x, y = columns(orbit)
    with_theme(theme_minimal(),
        markersize = 5) do
        
        X = x[1:300]
        Y = y[1:300]
        fig = Figure();
        ax_ts = Axis(fig[1, 1:2], xlabel = "Time (t)", ylabel = "Value")
        scatterlines!(ax_ts, (X .- mean(X)) ./ std(X), label = "x")
        scatterlines!(ax_ts, (Y .- mean(Y)) ./ std(Y), label = "y")
        axislegend()

        ax_hm = Axis(fig[2, 1:2], xlabel = "a", ylabel = "b")
        ax_hm.yticks = (1:length(as), string.([i % 2 == 0 ? as[i] : "" for i = 1:length(as)]))
        ax_hm.xticks = (1:length(bs), string.([i % 2 == 0 ? bs[i] : "" for i = 1:length(bs)]))
        hm = heatmap!(ax_hm, Δ,  colormap = :viridis)
        Colorbar(fig[2, 3], hm; label = "Δ' = ρ(ŷ | yx) - ρ(x̂ | xy)")
        fig
    end
end

reproduce_figure_8_mccraken()
```

We haven't used as many parameter combinations as [McCracken2014](@citet) did, 
but we get a figure that looks roughly similar to theirs.

As expected, ``\Delta < 0`` for all parameter combinations, implying that ``X`` "PAI drives" ``Y``.

### Spatial cross mapping

```@example example_spatial_cross_mapping
using CausalityTools
using DynamicalSystemsBase
using CairoMakie


function eom_tentmap(dx, x, p, n)
    x = x[1]
    μ = p[1]
    dx[1] = x < 0.5 ? μ*x : μ*(1 - x)

    return
end

function tentmap(u₀ = rand(); μ = 1.98)
    DiscreteDynamicalSystem(eom_tentmap, [u₀], [μ])
end 

npts = 2000
sys = tentmap(μ = 1.98)
x, tinds = trajectory(sys, npts , Ttr = 1000)
x = diff(first(columns(x)))

n = length(x)
k = 7 # number of time steps into the future to predict
nmax = n - (d*τ +1) - k
d, τ = 3, 1
function train_and_pred(k, n, d, τ)
    training = 1:(nmax ÷ 2)
    prediction = ((nmax ÷ 2) + 1):nmax
    return training, prediction
end
training, prediction = train_and_pred(1, n, d, τ)
x_1, x̃_1 = simplex_predictions(x, 1, d = d, τ = τ, training = training, prediction = prediction)

training, prediction = train_and_pred(7, n, d, τ)
x_7, x̃_7 = simplex_predictions(x, 7, d = d, τ = τ, training = training, prediction = prediction)

f = Figure(); ax = Axis(f[1, 1], xlabel = "Observed values", ylabel = "Predicted values")
scatter!(ax, x_1, x̃_1, label = "k = 1", marker = :circle)
scatter!(ax, x_7, x̃_7, label = "k = 7", marker = :circle)
axislegend()
f
```


There is high correlation between observed and predicted values when predicting only one time step (`k = 1`)
into the future. As `k` increases, the performance drops off. Let's investigate this systematically.

```@example example_spatial_cross_mapping
kmax = 20
cors = zeros(kmax)
for k = 1:kmax
    training, prediction = train_and_pred(k, n, d, τ)
    X, X̃ = simplex_predictions(x, k, d = d, τ = τ, training = training, prediction = prediction)
    cors[k] = cor(X, X̃)
end

f = Figure(); 
ax = Axis(f[1, 1]; 
    xlabel = "Prediction time (k)", ylabel = "Correlation coefficient (ρ)",
     limits = ((0.5, kmax + 0.5), (-1.1, 1.1))
)
#hlines!(ax, [0], linestyle = :dash, label = "", color = :grey)
scatter!(ax,  cors)
f
```


The correlation between observed and predicted values is near perfect until `k = 3`, and then rapidly 
drops off as `k` increases. At `k = 8`, there is virtually no correlation between observed and predicted values.
This means that, for this particular system, for this particular choice of embedding and choice of training/prediction sets, the predictability of the system is limited to about 4 or 5 time steps into the future (if you want good predictions). 

The main point of Sugihara & May's paper was that this drop-off of prediction accuracy with `k` is characteristic of chaotic systems, and can be used to distinguish chaos from regular behaviour in time series.

Let's demonstrate this by also investigating how the correlation between observed and predicted values behaves as a function of `k` for a regular, non-chaotic time series. We'll use a sine wave with additive noise.

```@example example_spatial_cross_mapping
using Distributions
𝒩 = Uniform(-0.5, 0.5)
xs = 0.0:1.0:2000.0
r = sin.(0.5 .* xs) .+ rand(𝒩, length(xs))
plot(r[1:200])

cors_sine = zeros(kmax)
for k = 1:kmax
    training, prediction = train_and_pred(k, n, d, τ)
    X, X̃ = simplex_predictions(r, k, d = d, τ = τ, training = training, prediction = prediction)
    cors_sine[k] = cor(X, X̃)
end


f = Figure(); 
ax = Axis(f[1, 1]; 
    xlabel = "Prediction time (k)", ylabel = "Correlation coefficient (ρ)",
     limits = ((0.5, kmax + 0.5), (-1.1, 1.1))
)
scatter!(ax, 1:kmax, cors, label = "tent map", marker = :star6, markersize = 10)
scatter!(ax, 1:kmax, cors_sine, label = "sine", marker = :hexagon, markersize = 8)
axislegend()
f
```

In contrast to the tent map, for which prediction accuracy drops off and stabilizes around zero for increasing `k`, the prediction accuracy is rather insensitive to the choice of `k` for the noisy sine time series. 


### Example: determining optimal embedding dimension

```@docs
delay_simplex
```

The simplex projection method can also be used to determine the optimal embedding dimension for a time series.
Given an embedding lag `τ`, we can embed a time series `x` for a range of embedding dimensions `d ∈ 2:dmax` and
compute the average prediction power over multiple `ks` using the simplex projection method.

Here, we compute the average prediction skills from `k=1` up to `k=10` time steps into the future, for 
embedding dimensions `d = 2:10`. We'll use a coupled Lorenz attractor system.

```@example example_spatial_cross_mapping
using CausalityTools
using CairoMakie
using DynamicalSystemsBase
using DelayEmbeddings

Base.@kwdef struct LorenzBidir6{V, CXY, CYX, A1, A2, A3, B1, B2, B3}
    xi::V = [0.1, 0.05, 0.2, 0.2, 0.25, 0.3]
    c_xy::CXY = 0.2
    c_yx::CYX = 0.2
    a₁::A1 = 10
    a₂::A2 = 28
    a₃::A3 = 8/3
    b₁::B1 = 10
    b₂::B2 = 28
    b₃::B3 = 9/3
end

function system(definition::LorenzBidir6)
    return ContinuousDynamicalSystem(eom_lorenzlorenzbidir6, definition.xi, definition)
end

@inline @inbounds function eom_lorenzlorenzbidir6(u, p, t)
    (; xi, c_xy, c_yx, a₁, a₂, a₃, b₁, b₂, b₃) = p
    x1, x2, x3, y1, y2, y3 = u

    dx1 = -a₁*(x1 - x2) + c_yx*(y1 - x1)
    dx2 = -x1*x3 + a₂*x1 - x2
    dx3 = x1*x2 - a₃*x3
    dy1 = -b₁*(y1 - y2) + c_xy*(x1 - y1)
    dy2 = -y1*y3 + b₂*y1 - y2
    dy3 = y1*y2 - b₃*y3

    return SVector{6}(dx1, dx2, dx3, dy1, dy2, dy3)
end

sys = system(LorenzBidir6())
T, Δt = 150, 0.05
lorenz, ts = trajectory(sys, T, Δt = Δt, Ttr = 100)
x1, x2, x3 = columns(lorenz)[1:3]

# Determine the optimal embedding delay
τ = estimate_delay(x1, "ac_zero")

# Compute average prediction skill for a range of dimensions. We'll average 
# the prediction skill over time steps `k = 1:10` for each dimension.
ds, ks = 2:10, 1:10
ρs = delay_simplex(x1, τ, ds = ds, ks = ks)

f = Figure(); ax = Axis(f[1,1]; xlabel = "Embedding dimension", ylabel = "ρ̄(observed, predicted")
scatterlines!(ax, ds, ρs, label = "", color = :black, marker = :star)
f
```

Based on the predictability criterion, the optimal embedding dimension, for this particular realization
of the first variable of the Lorenz system, seems to be 2.

## S-map

```@docs
smap
```

The s-map, or sequential locally weighted global map, was introduced in Sugihara (1994)[^Sugihara1994]. The s-map approximates the dynamics of a system as a locally weighted global map, with a tuning parameter ``\theta`` that controls the degree of nonlinearity in the model. For ``\theta = 0``, the model is the maximum likelihood global linear solution (of eq. 2 in Sugihara, 1994), and for increasing ``\theta > 0``, the model becomes increasingly nonlinear and localized (Sugihara, 1996)[^Sugihara1996].

When such a model has been constructed, it be used as prediction tool for out-of-sample points, and can be used to characterize nonlinearity in a time series (Sugihara, 1994).  Let's demonstrate with an example.

### Example: prediction power for the Lorenz system

In our implementation of `smap`, the input is a multivariate dataset - which can be a `StateSpaceSet` of either 
the raw variables of a multivariate dynamical system, or a `Dataset` containing an embedding of a single time series. 
Here, we'll show an example of the former.

Let's generate an example orbit from a the bidirectionally coupled set of Lorenz systems defined above. 
We'll select the first three variables for analysis.

```@example example_spatial_cross_mapping
using CausalityTools
using CairoMakie
using DynamicalSystemsBase
using Statistics


sys = system(LorenzBidir6())
T, Δt = 150, 0.05
lorenz, ts = trajectory(sys, T, Δt = Δt, Ttr = 100)
lorenz = lorenz[:, 1:3]
x1, x2, x3 = columns(lorenz)

f = Figure(); ax = Axis3(f[1,1]; xlabel = "x", ylabel = "y", zlabel = "z")
scatterlines!(ax, x1, x2, x3, marker = :circle, label = "", markersize = 2, alpha = 0.5, linewidth = 1)
f
```

Now, we compute the `k`-step forward predictions for `k` ranging from `1` to `15`. The tuning parameter `θ` 
varies from `0.0` (linear model) to `2.0` (strongly nonlinear model). Our goal is to see which model 
yields the best predictions across multiple `k`.

We'll use the first 500 points of the orbit to train the model. Then, using that model, we try to
predict the next 500 points (which are not part of the training set). 
Finally, we compute the correlation between the predicted values and the observed values, which measures
the model prediction skill. This procedure is repeated for each combination of `k` and `θ`.

```@example smap_lorenz
using CausalityTools
using CairoMakie
using DynamicalSystemsBase
using Statistics

function train_and_pred_smap(nmax)
    training = 1:(nmax ÷ 2)
    prediction = ((nmax ÷ 2) + 1):nmax
    return training, prediction
end

T, Δt = 30, 0.05
lorenz, ts = trajectory(sys, T, Δt = Δt, Ttr = 100)
ks, θs = 1:10, 0.0:0.5:2.0
n = length(lorenz)

# Compute correlations between predicted values `preds` and actual values `truths` 
# for all parameter combinations
cors = zeros(length(ks), length(θs))
for (i, k) in enumerate(ks)
    println("k=$k")
    for (j, θ) in enumerate(θs)
        println("θ=$θ")
        training, prediction = train_and_pred_smap(n)
        preds, truths = smap(lorenz, θ = θ, k = k, trainees = training, predictees = prediction)
        cors[i, j] = cor(preds, truths)
    end
end
cors

f = Figure(); 
ax = Axis(f[1,1]; 
    xlabel = "Prediction time (k)", ylabel = "cor(observed, predicted)",
    limits = (nothing, (-1.1, 1.1)),
)
markers = [:star :circle :square :star5 :hexagon :circle :star]
cols = [:black, :red, :blue, :green, :purple, :grey, :black]

labels = ["θ = $θ" for θ in θs]
for i = 1:length(θs)
    scatterlines!(ks, cors[:, i], marker = markers[i], color = cols[i], markersize = 5, label = labels[i],
        linewidth = i == 1 ? 5 : 2)
end
axislegend()
f
```

The nonlinear models (colored lines and symbols) far outperform the linear model (black line + stars).

Because the predictions for our system improves with increasingly nonlinear models, it indicates that our 
system has some inherent nonlinearity. This is, of course, correct, since our Lorenz system is chaotic.

A formal way to test the presence of nonlinearity is, for example, to define the null hypothesis 
"H0: predictions do not improve when using an equivalent nonlinear versus a linear model" (equivalent in the sense  that the only parameter is `θ`) or, equivalently, "`H0: ρ_linear = ρ_nonlinear`". If predictions do in fact improve, we
instead accept the alternative hypothesis that prediction *do* improve when using nonlinear models versus using linear models. This can be formally tested using a z-statistic [^Sugihara1994].

[^Sugihara1994]: Sugihara, G. (1994). Nonlinear forecasting for the classification of natural time series. Philosophical Transactions of the Royal Society of London. Series A: Physical and Engineering Sciences, 348(1688), 477-495.
[^Sugihara1996]: Sugihara, George, et al. "Nonlinear control of heart rate variability in human infants." Proceedings of the National Academy of Sciences 93.6 (1996): 2608-2613.