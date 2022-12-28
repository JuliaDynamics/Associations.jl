
## Continuous/differential entropy

Continuous (differential) entropies are defined by an *integral*, and are 
related to but don't share all the same properties as discrete entropies.
For example, Shannon differential entropy may even be *negative* for 
some distributions.

Continuous entropies must be estimated using some form of "plug-in" estimator. For example, the Shannon differential entropy for a random variable $X$ with support $\mathcal{X}$ is defined as

$$h(x) = \mathbb{E}[-\log{(f(X))}] = -\int_{\mathcal{X}}f(x) \log f(x) dx.$$

There are several ways of estimating this integral from observed data,
using what is called "plug-in" estimators. A common plug-in estimator is the resubstitution estimator

$$\hat{H}(x) = -\frac{1}{N}\sum_{i=1}^N \log{(\hat{p}(X_i))},$$

where $\hat{p}$ is estimated using the samples $X_1, X_2, \ldots, X_N$, is a plug-in estimator for Shannon differential entropy.

Subtypes of [`DifferentialEntropyEstimator`](@ref)s use various forms of plug-in estimators to estimate differential entropy. For example, [`Kraskov`](@ref) estimates *Shannon* differential entropy. [`LeonenkoProzantoSavani`](@ref), on the other hand, estimates both [`Shannon`](@ref), [`Renyi`](@ref) and 
[`Tsallis`](@ref) differential entropy.

!!! note "Plug-in estimators for differential entropy"
    When using [`entropy`](@ref) with a [`ProbabilityEstimator`](@ref), it is always the discrete entropy that is computed. When using [`entropy`](@ref) with an [`DifferentialEntropyEstimator`](@ref), it is the *differential* entropy that is computed.


## Generalized entropies

There exists a multitude of various entropy measures in the scientific literature, which
appear either in discrete form, differential/continuous form, or both.

* Discrete entropies are simply functions of sums over
[probability mass functions (pmf)](https://en.wikipedia.org/wiki/Probability_mass_function).
Hence, every [`ProbabilitiesEstimator`](@ref) yields a naive plug-in estimator for
generalized entropy (i.e. just plug the probabilities into the relevant entropy formulas).
No bias correction is currently performed for discrete estimators.

* Continuous (differential) entropies are functions of
[probability density functions (pdf)](https://en.wikipedia.org/wiki/Probability_density_function).
Therefore, continuous entropy estimators approximate *integrals* instead of sums, as
in the discrete case, and boils down to density function estimation. Every
[`DifferentialEntropyEstimator`](@ref)s has a slightly different way of estimating densities, 
hence yielding slightly different differential entropy estimates.
