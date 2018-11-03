# What is a surrogate?
To generate a surrogate for a given data set means generating another data series that preserves some statistical or dynamical property of the original data set.

Consider the following example, where we generate different types of surrogate realizations for a dataset of sunspot occurrences.

## Example: surrogate realizations of sunspot data

```@example
# Load the sunspot data which was downloaded from
# http://www.sidc.be/silso/datafiles#total on September 9th 2018
sunspots = readdlm("../../data/SN_d_tot_V2.0.txt")

# Each row in this array represents a single day. The fifth column in the
# number of sunspots. To generate surrogates, we'll use the sunspot number
# from August 31st 2017 to August 31st 2018.
sunspots_2017 = sunspots[end-365:end, :]

using Plots
using TimeseriesSurrogates

lf = Plots.font("Helvetica", 10)

p1 = plot(sunspots_2017[:, 5], lc = :black,
        title = "Sunspot data", titlefont = lf)
p2 = plot(randomshuffle(float.(sunspots_2017[:, 5])),
        title = "Random shuffle surrogate", titlefont = lf)
p3 = plot(randomphases(float.(sunspots_2017[:, 5])),
        title = "Random phase Fourier surrogate", titlefont = lf)
p4 = plot(randomamplitudes(float.(sunspots_2017[:, 5])),
        title = "Random amplitude Fourier surrogate", titlefont = lf)
p5 = plot(iaaft(float.(sunspots_2017[:, 5])),
        title = "IAAFT surrogate", titlefont = lf)
p6 = plot(aaft(float.(sunspots_2017[:, 5])),
        title = "AAFT surrogate", titlefont = lf)
plot(p1, p2, p3, p4, p5, p6, layout = (3, 2),
        size = (700, 600), legend = false)
savefig("sunspot_surrogates.svg"); nothing #hide  
```
![](sunspot_surrogates.svg)

## Types of surrogate methods
There are two types of surrogate realizations of a data set:

1. **Constrained realizations**: Realizations that shuffle around the original values of the data, either destroying or retaining certain statistical properties of the original data. One example is the random shuffle surrogate, where the values of
a time series is shuffled randomly, preserving the amplitude distribution of the original time series, but destroys any other dynamical information ([`randomshuffle`](@ref)). Other examples: amplitude-adjusted Fourier transform surrogates ([`aaft`](@ref)), and iterated amplitude-adjusted Fourier transform  surrogates ([`iaaft`](@ref)).
2. **Typical realizations**: Realizations created by fitting some model to the data, then generating the surrogate realizations from that model. Examples: [`randomamplitudes`](@ref) and [`randomphases`](@ref).

In this package, surrogates are used, among other things, for null-hypothesis testing when applying the causality detection algorithms.

The surrogate functions here are re-exported from [`TimeseriesSurrogates.jl documentation`](https://github.com/kahaaga/TimeseriesSurrogates.jl). Please consult the [`TimeseriesSurrogates.jl documentation`](https://kahaaga.github.io/TimeseriesSurrogates.jl/latest/) for the complete API.
