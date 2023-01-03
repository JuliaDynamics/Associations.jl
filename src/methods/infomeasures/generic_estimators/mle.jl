"""
    MaximumLikelihoodEstimator <: Estimator
    MaximumLikelihoodEstimator()
    MLE() # alias

A maximum likelihood (ML), "plug-in" or "naive" estimator.

## Description

In the context of information theory, a maximum likelihood estimator
simply plugs in the maximum likelihood estimate of some probabilities
(i.e. the relative frequencies of occurrences of the observed outcomes;
Arora et al., 2022)) into some function of probabilities.
"""
struct MaximumLikelihoodEstimator <: Estimator end
