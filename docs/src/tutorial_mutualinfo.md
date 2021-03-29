# Mutual information on synthetic systems

In this example we generate realizations of two different systems where we know the strength of coupling between the variables. Our aim is to compute mutual information $I(X; Y)$ between time series of each variable and assess how the magnitude of $I(X; Y)$ changes as we change the strength of coupling between $X$ and $Y$.

## Defining the systems

Here we implement two of the example systems that come with the [CausalityTools.jl](https://github.com/JuliaDynamics/CausalityTools.jl) package:

* A stochastic system consisting of two unidirectionally coupled first-order autoregressive processes ([`ar1_unidir`](@ref))
* A deterministic, chaotic system consisting of two unidirectionally coupled logistic maps ([`logistic2_unidir`](@ref))

We use the default input parameter values (see [`ar1_unidir`](@ref) and [`logistic2_unidir`](@ref) for details) and below we toggle only the random initial conditions and the coupling strength parameter `c_xy`. For each value of `c_xy` we generate 1,000 unique realizations of the system and obtain 500-point time series of the coupled variables.

## Estimating mutual information

Here we use the binning-based [`VisitationFrequency`](@ref) estimator. We summarize the distribution of $I(X; Y)$ values across all realizations using the median and quantiles encompassing 95\% of the values.

```@example
using CausalityTools, Statistics, Plots

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
mi_estimator = VisitationFrequency(b)

for i in 1 : length(c)
    
    tmp = zeros(n_realizations, 2)
    
    for k in 1 : n_realizations
        
        # Obtain time series realizations of the two 2D systems 
        # for a given coupling strength and random initial conditions
        lmap = trajectory(logistic2_unidir(u₀ = rand(2), c_xy = c[i]), npts - 1, Ttr = 1000)
        ar1 = trajectory(ar1_unidir(u₀ = rand(2), c_xy = c[i]), npts - 1)
        
        # Compute the MI between the two coupled components of each system
        tmp[k, 1] = mutualinfo(lmap[:, 1], lmap[:, 2], mi_estimator)
        tmp[k, 2] = mutualinfo(ar1[:, 1], ar1[:, 2], mi_estimator)
    end
    
    # Compute lower, middle, and upper quantiles of MI for each coupling strength
    mi[i, :, 1] = quantile(tmp[:, 1], [0.025, 0.5, 0.975])
    mi[i, :, 2] = quantile(tmp[:, 2], [0.025, 0.5, 0.975])
end

# Plot distribution of MI values as a function of coupling strength for both systems
Plots.plot(c, mi[:, 2, 1], label = "2D chaotic logistic maps", lc = "black",
    ribbon = (mi[:, 2, 1] - mi[:, 1, 1], mi[:, 3, 1] - mi[:, 2, 1]), c = "black", fillalpha = 0.3,
    legend = :topleft)
Plots.plot!(c, mi[:, 2, 2], label = "2D order-1 autoregressive", lc = "red",
    ribbon = (mi[:, 2, 2] - mi[:, 1, 2], mi[:, 3, 2] - mi[:, 2, 2]), c = "red", fillalpha = 0.3)
Plots.xlabel!("Coupling strength")
Plots.ylabel!("Mutual information")

Plots.savefig("./mutualinfo_example.png")
```

![png](./mutualinfo_example.png)

As expected, $I(X; Y)$ increases with coupling strength in a system-specific manner.
