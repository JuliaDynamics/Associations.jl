# [`PairwiseAsymmetricInference`](@ref)

## Reproducing McCracken & Weigel (2014)

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
