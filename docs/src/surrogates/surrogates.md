# What is a surrogate?
To generate a surrogate for a given data set means generating another data series that preserves some statistical or dynamical property of the original data set.

One example is a random shuffle surrogate, where the values of
a time series is shuffled randomly. This preserves the amplitude distribution of the original time series, but destroys any other dynamical information (for example, static correlations or autocorrelation structure).

In addition to the random shuffle surrogates, the `TimeseriesSurrogates.jl` package provides Fourier phase surrogates, Fourier amplitude surrogates, amplitude-adjusted Fourier transform (AAFT), and iterated amplitude-adjusted Fourier transform (IAAFT) surrogates. These are all re-exported to `CausalityTools.jl` and can be used for null-hypothesis testing when applying
any of the causality detection algorithms.

Examples of surrogate generation can be found in the menu. Please consult the [`TimeseriesSurrogates.jl documentation`](https://kahaaga.github.io/TimeseriesSurrogates.jl/latest/) for the complete API.
