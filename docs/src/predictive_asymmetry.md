# Predictive asymmetry

```@docs
predictive_asymmetry
```

## Example

Here, we'll compute the predictive asymmetry on 100 different sets of random time series. Because there is no
dynamical coupling between the time series, we expect the predictive asymmetry to be zero.

We'll use a visitation frequency estimator.

```@example predictiveasymmetry
using CausalityTools, Plots, StatsBase

# Define prediction lags and estimator
ηs = 1:10
est = VisitationFrequency(RectangularBinning(5)) # guess that 5 bins along each coordinate axis is sufficient

nreps = 100
𝔸xy = zeros(nreps, length(ηs))
𝔸yx = zeros(nreps, length(ηs))

for i = 1:nreps
    # Some example time series
    x, y = rand(1000), rand(1000)
    𝔸xy[i, :] = predictive_asymmetry(x, y, est, ηs) 
    𝔸yx[i, :] = predictive_asymmetry(y, x, est, ηs) 
end

plot()
plot!(ηs, dropdims(mean(𝔸xy, dims = 1), dims = 1), label = "X to Y", c = :black)
plot!(ηs, dropdims(mean(𝔸yx, dims = 1), dims = 1), label = "Y to X", c = :red)
ylims!((-0.1, 0.1))
xlabel!("Prediction lag (η)")
ylabel!("Predictive asymmetry (bits)")
savefig("random-pa.svg"); nothing # hide
```

![](random-pa.svg)

As expected, because there is no dynamical coupling between the variables, the predictive asymmetry is around zero for all prediction lags.