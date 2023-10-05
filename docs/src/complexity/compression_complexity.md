# [Compression complexity](@id compression_complexity)

## Interface

```@docs
compression_complexity
```

## Algorithms

```@docs
EffortToCompress
```

## Example

Below is a reproduction of figure 3 in Nagaraj et al. (2013)[^Nagaraj2013], which compares the (approximate) Lyapunov exponent and compression complexity (ETC) of 200-pt long logistic map time series with varying bifurcation parameters.

For ETC, the time series is symbolized to a binary time series before analysis. Lyapunov exponents are estimated directly on the (non-symbolized) time series.

```@example
using CausalityTools, Plots, StatsBase, Measures

# Approximation to the Lyapunov exponent for a time series x, from Nagaraj (2013)
function lyapunov(x, a)
    L = length(x)
    lyp = 0.0
    for i in 2:length(x)
        lyp += log(2, abs(a - 2*a*x[i]))
    end
    
    return lyp / L
end


coeffs = 3.5:0.0001:4.0
ls = zeros(length(coeffs))
etcs = zeros(length(coeffs))
for (i, a) in enumerate(coeffs)
    sys = logistic2_unidir(r₁ = a, c_xy = 0.0, σ = 0.0)
    # Generate time series and compute approximate Lyapunov exponents
    x = trajectory(sys, 200, Ttr = 10000)[:, 1]
    ls[i] = lyapunov(x, a)

    # Symbolize and compute effort-to-compress
    y = [xᵢ > 0.5 ? 1 : 0 for xᵢ in x] 
    etcs[i] = compression_complexity(y, EffortToCompress(normalize = false))
end

plot(xlabel = "a", 
    legend = :topleft, 
    right_margin = 10mm, 
    bottom_margin = 5mm)
plot!(coeffs, ls .- mean(ls), label = "Lyapunov exponent", c = :red)
plot!(twinx(), coeffs, etcs, label = "ETC", axis = :right, c = :black)
savefig("compression_complexity_etc.svg") # hide
```

![](compression_complexity_etc.svg)

[^Nagaraj2013]: Nagaraj, N., Balasubramanian, K., & Dey, S. (2013). A new complexity measure for time series analysis and classification. The European Physical Journal Special Topics, 222(3), 847-860.
