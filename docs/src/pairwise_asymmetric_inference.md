# Pairwise asymmetric inference

```@docs
pai
```

## Example: nonlinear system

Let's try to reproduce figure 8 in McCracken & Weigel (2014)[^McCracken2014]. We'll start by defining the their example B (equations 6-7). This system consists of two
variables ``X`` and ``Y``, where ``X`` drives ``Y``.

```@example pai_ex
using CausalityTools, DynamicalSystems, Plots, StatsBase, Statistics, Distributions; gr()

function eom_nonlinear_sindriver(dx, x, p, n)
    a, b, c, t, Δt = (p...,)
    x, y = x[1], x[2]
    𝒩 = Normal(0, 1)
    
    dx[1] = sin(t)
    dx[2] = a*x * (1 - b*x) + c*rand(𝒩)
    p[end-1] += 1 # update t

    return
end

function nonlinear_sindriver(;u₀ = rand(2), a = 1.0, b = 1.0, c = 2.0, Δt = 1)
    DiscreteDynamicalSystem(eom_nonlinear_sindriver, u₀, [a, b, c, 0, Δt])
end

# Create a system of nonidentical logistic maps where coupling from variable x to variable y
# is stronger than vice versa.
sys = nonlinear_sindriver(a = 1.0, b = 1.0, c = 2.0)
npts = 100
orbit = trajectory(sys, npts, Ttr = 10000)
x, y = columns(orbit);
plot(xlabel = "Time step", ylabel = "Value")
# Standardize and plot data
plot!((x .- mean(x)) ./ std(x), label = "x")
plot!((y .- mean(y)) ./ std(y), label = "y")
savefig("pai_ts.svg") # hide
```

![](pai_ts.svg)

Now, let's generate such time series for many different values of the parameters `a` and `b`, and compute PAI for fixed `p = 2.0`. This will replicate the upper right panel of figure 8 in the original paper.

```@example pai_ex
as = 0.25:0.25:4.0
bs = 0.25:0.25:4.0

pai_xys = zeros(length(as), length(bs))
pai_yxs = zeros(length(as), length(bs))
c = 2.0
npts = 2000
d, τ = 2, 1
for (i, a) in enumerate(as)
    for (j, b) in enumerate(bs)
        s = nonlinear_sindriver(a = a, b = a, c = c)
        orbit = trajectory(s, npts, Ttr = 10000)
        X, Y = columns(orbit)
        # Use the segment bootstrap estimator, take the mean of 50 reps over segments of # length L = 200
        pai_xys[i, j] = pai(X, Y, d, τ, :segment, L = 200, nreps = 50) |> mean
        pai_yxs[i, j] = pai(Y, X, d, τ, :segment, L = 200, nreps = 50) |> mean
    end
end
```

Now that we have computed the PAI in both directions, we define a measure of directionality as the difference between PAI in the ``X`` to ``Y`` direction and in the ``Y`` to ``X`` direction, so that if ``X`` drives ``Y``, then ``\\Delta < 0``.

```@example pai_ex
Δ = pai_xys .- pai_yxs

clr = cgrad(:magma, categorical = true)
plot(xlabel = "a", ylabel = "b")
yticks!((1:length(as), string.(as)))
xticks!((1:length(bs), string.(bs)))
heatmap!(Δ, c = clr, logscale = true)
savefig("heatmap_pai.svg") # hide
```

![](heatmap_pai.svg)

As expected, ``\Delta < 0`` for all parameter combinations, implying that ``X`` "PAI drives" ``Y``.