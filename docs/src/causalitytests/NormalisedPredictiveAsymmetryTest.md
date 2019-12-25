# [`NormalisedPredictiveAsymmetryTest`](@id normalised_predictive_asymmetry)

```@docs
NormalisedPredictiveAsymmetryTest
```

## [Example](@id example_NormalisedPredictiveAsymmetryTest)

### Preparations

For this example, we'll use the built-in `logistic4` system, which consists of 
four unidirectionally coupled logistic maps coupled ``x \to y \to z \to w``.
With default values, the system blows up for some initial conditions, so we'll 
start by just making a thin wrapper that finds us a good orbit.

```@example NormalisedPredictiveAsymmetryTest_logistic4
using CausalityTools, DynamicalSystems
using Plots, LaTeXStrings; gr()

"""
    good_orbit(f::Function, npts::Int; Ttr::Int = 100, max_attempts = 1000)

Draw a good orbit from a dynamical system defined by `f` when called. Useful
when some initial conditions yield attractors while others do not. CAREFUL:
this function is defined recursively.
"""
function good_orbit(f::Function, npts::Int; Ttr::Int = 100, max_attempts = 1000)
    sys = f()
    ts = columns(trajectory(sys, npts, Ttr = Ttr))
    
    tries = 0
    while tries < max_attempts
        tries += 1
        if all(isfinite.(hcat(ts...)))
            return ts
        end
    end
    println("Did not manage to find good orbits within $(max_attempts) attempts. Re-initializing.")
    
    good_orbit(f, npts, Ttr = Ttr, max_attempts = max_attempts)
end
```

### Generate time series

Okay, now we're ready to generate some time series.

```@example NormalisedPredictiveAsymmetryTest_logistic4
# Example data
sys = logistic4
npts = 400

# Some initial conditons blow up, so iterate until we find a good one.
x, y, z, w = good_orbit(sys, npts, Ttr = 300)

p_ts = plot(xlabel = "Time step", ylabel = "Value", size = (382*2, 300),
    guidefont = font(13), tickfont = font(13), legendfont = font(13))
plot!(x, label = "x")
plot!(y, label = "y")
plot!(z, label = "z")
plot!(w, label = "w")
```

### Test setup

Now, define a normalised predictive asymmetry test using a visitation frequency test 
to get our predictions for lags `[-10:1; 1:10]`, symmetrically around the zero lag. 
To construct the  [partition](@ref partition_rectangular) of our 
[state space reconstructions](@ref custom_delay_reconstruction), we'll use the 
heuristic from [^1] that the number of intervals along each coordinate axis 
should not exceed ``N^{\left(\frac{1}{D+1}\right)}``, where ``N`` is the number of points in 
the time series and ``D`` is the total dimension of the generalised reconstruction
used for the transfer entropy analyses.

```@example NormalisedPredictiveAsymmetryTest_logistic4
# Test setup
η_max = 10
k = 1; l = 1; m = 1 # embedding parameters, total dim = k + l + m
ηs = [-η_max:-1; 1:η_max]
bin = RectangularBinning(floor(Int, npts^(1/(k + l + m + 1))))
test = VisitationFrequencyTest(ηs = ηs, binning = bin)
pa_test = NormalisedPredictiveAsymmetryTest(test, f = 1.0)
```

That's all we need. Now, let's compute the predictive asymmetry statistic and 
see whether we find the expected causalities.

```@example NormalisedPredictiveAsymmetryTest_logistic4
pas_xy = causality(x, y, pa_test)
pas_yx = causality(y, x, pa_test)
pas_yz = causality(y, z, pa_test)
pas_zy = causality(z, y, pa_test)
pas_zw = causality(z, w, pa_test)
pas_wz = causality(w, z, pa_test)
asymmetries = [pas_xy pas_yx pas_yz pas_zy pas_zw pas_wz]
```

### Plot time series and results

Let's plot the results.

```@example NormalisedPredictiveAsymmetryTest_logistic4
ymax = maximum(abs.(asymmetries)) * 1.05
p_pa = plot(xlabel = "\$ \\eta \$", ylabel = "\$ \\mathcal{A}(\\eta) \$",
    ylims = (-ymax, ymax),  size = (382*2, 550), legend = :topleft,
    fg_legend = :transparent, bg_legend = :transparent,
    guidefont = font(13), tickfont = font(13), legendfont = font(12),
    xlims = (0.9, η_max+0.1))
hline!([1], ls = :dash, lc = :grey, label = "")
plot!(ηs[ηs .> 0], pas_xy, c = :black, label = "\$ x \\to y \$")
plot!(ηs[ηs .> 0], pas_yz, c = :red, label = "\$ y \\to z \$")
plot!(ηs[ηs .> 0], pas_zw, c = :blue, label = "\$ z \\to w \$")
plot!(ηs[ηs .> 0], pas_yx, c = :black, ls = :dash, label = "\$ y \\to x \$")
plot!(ηs[ηs .> 0], pas_zy, c = :red, ls = :dash, label = "\$ z \\to y \$")
plot!(ηs[ηs .> 0], pas_wz, c = :blue, ls = :dash, label = "\$ w \\to z \$")
```

### Discussion

If the predictive asymmetry picks up on the causal interactions, it should be positive 
for the directions we know to be causal and negative in the directions that are known to 
be non-causal. Running this example for multiple initial conditions of the system, we 
find that this is very consistently the case!

[^1]: 
    Krakovská, A., Jakubík, J., Chvosteková, M., Coufal, D., Jajcay, N., & Paluš, M. (2018). Comparison of six methods for the detection of causality in a bivariate time series. Physical Review E, 97(4), 042207.
    [https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.042207](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.042207)