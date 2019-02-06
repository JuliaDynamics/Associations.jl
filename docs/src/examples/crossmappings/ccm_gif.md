## Disclaimer:
This function is not an integrated part of the `CausalityTools.jl` package
and therefore neither formally documented, nor covered by tests or input validation.
Use it at your own risk!

## Description

`make_ccm_gif` uses `convergentcrossmap` to create plots of two time series vectors `x`and `y` at increasing `ts_lengths` and compile them into a gif. The final plot is a composite of four subplots:

1. A driver plot, assuming `x` as the driver time series.
2. A response plot, assuming `y` as response time series.
3. A skill plot, showing prediction skill ``\rho`` for ``x \rightarrow y`` and ``y \rightarrow x``.
4. An embedding plot, showing the three-dimensional attractor for the response time series `y` using the embedding lag `τ`.

The number of generated plots depends on the length of `ts_lengths`; one plot is created per element.

`make_ccm_gif` accepts all `kwargs` taken by [`convergentcrossmap`](../../crossmappings/crossmapping.md). Other `kwargs` taken by `make_ccm_gif` specifically:
`fps::Int` sets the fps for the saved gif, defaults to `5`. `save_anim::Bool` saves the animation as a
`.jld2` file if `true`, but defaults to `false`. `τ::Int` sets embedding lag necessary to create
the 3D embedding plot, defaults to `1`. Once run, the function automatically saves output as `tmp.gif` in the
current working directory.

## Code

This is the code used to generate all animated ccm plots in the `CausalityTools.jl` documentation.

```julia
using PyPlot
using PyCall
using Plots
using Distributions
using CausalityTools
using Measures
using JLD2
using Dates

function composite_ccm_plot!(x, y, driver_plot, response_plot, embed_plot, skill_plot;
                           Lmin, Lmax, τ, kwargs...)

    # Calculate cross mapping skills
    average_xy, uncertainty_xy = convergentcrossmap(x, y, [Lmin, Lmax], dim = 3, τ = τ, kwargs...)
    average_yx, uncertainty_yx = convergentcrossmap(y, x, [Lmin, Lmax], dim = 3, τ = τ, kwargs...)

    # Driver and response plots
    Plots.plot!(driver_plot, Lmin:Lmax, x[Lmin:Lmax], lc = :red, label = "", lw = 0.3)
    Plots.plot!(response_plot, Lmin:Lmax, y[Lmin:Lmax], lc = :blue, label = "", lw = 0.3)

    # Skill plot
    Plots.plot!(skill_plot, [Lmin, Lmax], average_xy, label = "", lc = :red, lw = 0.6, ls = :solid)
    Plots.plot!(skill_plot, [Lmin, Lmax], uncertainty_xy, lc = :red, ls = :dash, label = "", lw = 0.6)
    Plots.plot!(skill_plot, [Lmin, Lmax], average_yx, lc = :blue, ls = :solid, label = "", lw = 0.6)
    Plots.plot!(skill_plot, [Lmin, Lmax], uncertainty_yx, lc = :blue, ls = :dash, label = "", lw = 0.6)

    # 3D embedding
    embedding = StateSpaceReconstruction.embed([y], [1, 1, 1], [0, -1, -2]).points
    Plots.scatter3d!(embed_plot, ([embedding[i, Lmin:Lmax] for i in 1:3]...,),
                label = "", mc = :black,  ms = 0.6, mα = 0.5)
end


function make_ccm_gif(x, y, ts_lengths; fps::Int = 5, save_anim::Bool = false, τ::Int = 1, kwargs...)

    lag1 = τ; lag2 = 2*lag1
    x_extremes = (minimum(x), maximum(x)) # x extreme values
    y_extremes = (minimum(y), maximum(y)) # y extreme values
    pyplot()

    # Define plots
    driver_plot = Plots.plot(x[1:minimum(ts_lengths)], ylabel = "Driver", lc = :red, label = "",
                                    xaxis = false, lw = 0.3,
                                    xlims = (0, maximum(ts_lengths)), ylims = x_extremes)
    response_plot = Plots.plot(y[1:minimum(ts_lengths)], ylabel = "Response", xlabel = "Time series length (L)",
                                    lc = :blue, label = "", lw = 0.3,
                                    xlims = (0, maximum(ts_lengths)), ylims = y_extremes)
    skill_plot = Plots.plot(ylabel = "ρ", xlabel = "Time series length (L)",
                                lc = :red, label = "", lw = 0.6,
                                xlims = (0, maximum(ts_lengths)), ylims = (0, 1))
    embed_plot = Plots.plot(label = "", mc = :blue, ms = 0.6, mα = 0.5,
                                xlabel = "response(t)", ylabel = "response(t-$lag1)", zlabel = "response(t-$lag2)",
                                right_margin = 6mm, tickfont = font(7),
                                xlims = y_extremes, ylims = y_extremes, zlims = y_extremes)
    emptyplot = Plots.plot(legend = false, grid = false, foreground_color_subplot = :white, mα = 0)

    # Create animation
    anime = @animate for (i, L) in enumerate(ts_lengths[2:end])
        if i == 1
            Plots.plot!(skill_plot, [-2, -1], [-1, -0.5], label = "x → y", lc = :red, lw = 0.6, ls = :solid)
            Plots.plot!(skill_plot, [-2, -1], [-1, -0.5], label = "y → x", lc = :blue, lw = 0.6, ls = :solid)
        end

        composite_ccm_plot!(x, y, driver_plot, response_plot,
                                    embed_plot, skill_plot,
                                    Lmin = ts_lengths[i], Lmax = L, τ = τ; kwargs...)
        driver_response_plot = Plots.plot(driver_plot, response_plot, layout = (2, 1))
        skill_attractor_plot = Plots.plot(skill_plot, embed_plot,
                                            emptyplot, layout = @layout [a{0.47w} b{0.51w} d{0.02w}])
        Plots.plot(driver_response_plot, skill_attractor_plot, layout = (2, 1), guidefont = font(8))
    end

    if save_anim
        JLD2.@save string("ccm_anim_", Dates.now(Dates.UTC), ".jld2") anime
    end

    # Return
    gif(anime, fps = fps)

end
```

Sample function call:

```julia
system = CausalityTools.Systems.ar1()
tra = trajectory(system, 1000-1)
x, y = tra[:, 1], tra[:, 2]
ts_lengths = 50:10:500

make_ccm_gif(x, y, ts_lengths)
```
