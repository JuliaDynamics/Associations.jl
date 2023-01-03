"""
    Jackknife <: Estimator
    Jackknife()
    Jackknife() # alias

A jackknife estimator.

## Description

A jackknife estimate of some quantity ``F`` over observations
``\\bar{X}_{1:n} = \\{x_1, x_2, \\ldots, x_n \\}``
is the average estimate, after estimating ``F`` over all subsets of
... # Todo
"""
struct MaximumLikelihoodEstimator <: Estimator end
