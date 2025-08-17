# Associations.jl

```@docs
Associations
```

## Latest news

!!! info "Package rename"
    The package has been *renamed* from CausalityTools.jl to Associations.jl.

Associations.jl has been updated to v4! 

This update includes a number of breaking changes, several of which are *not* backwards compatible.
These are done to ensure compatibility with 
[ComplexityMeasures.jl v3](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/), which provides discretization functionality that we use here.

Important changes are:
- Convenience methods have been removed completely. Use [`association`](@ref) instead.
- Example systems have been removed.
- The syntax for computing an association has changed. Estimators now always *contain the definition it estimates*. For example, `association(MIShannon(), KSG1(), x, y)` is now `association(KSG1(MIShannon()), x, y)`.
- `SurrogateTest` has been renamed to [`SurrogateAssociationTest`](@ref). 
- See the CHANGELOG.md for a complete list of changes.

## Getting started

The quickest way to get going with the package is to check out the examples in the left-hand menu.

!!! info
    To make it easier to navigate the extensive documentation, all documentation strings are 
    collapsed by default. Click the arrow icon in 
    the top toolbar to expand/collapse the docstrings in a page.


## Documentation content 

- [Association measures](@ref association_measures) lists all implemented association measures and their estimators.
- [Independence testing](@ref independence_testing) lists all implemented ways of determining if an association between datasets is "significant".
- [Causal inference](@ref causal_graphs) lists all methods of inferring association networks
  (also called "network graphs" and "causal graphs") between multiple variables.
- Numerous examples for [association measure estimation](@ref examples_associations), 
  [independence testing](@ref examples_independence), and 
  [network inference](@ref examples_network_inference).

## Input data for Associations.jl

Input data for Associations.jl are given as:

- Univariate *timeseries*, which are given as standard Julia `Vector`s.
- Multivariate timeseries, *StateSpaceSets*, or *state space sets*, which are given as
    [`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet)s. Many methods convert *timeseries* inputs to [`StateSpaceSet`](@extref StateSpaceSets.StateSpaceSet)
    for faster internal computations.
- Categorical data can be used with [`JointProbabilities`](@ref) to compute various
    information theoretic measures and is represented using any iterable whose elements
    can be any arbitrarily complex data type (as long as it's hashable), for example
    `Vector{String}`, `{Vector{Int}}`, or `Vector{Tuple{Int, String}}`.

## Maintainers and contributors

The Associations.jl software is maintained by
[Kristian Agasøster Haaga](https://github.com/kahaaga), who also curates and writes this
documentation. Significant contributions to the API and documentation design has been
made by [George Datseris](https://github.com/Datseris), which also co-authors
[ComplexityMeasures.jl](https://github.com/JuliaDynamics/ComplexityMeasures.jl), which
we develop in tandem with this package.

A complete list of contributors to this repo are listed on the main Github page. Some
important contributions are:

- [Norbert Genera](https://github.com/norbertgerena) contributed bug reports and
    investigations that led to subsequent improvements for the pairwise asymmetric
    inference algorithm and an improved cross mapping API.
- [David Diego](https://www.researchgate.net/profile/David-Diego)'s contributions were
    invaluable in the initial stages of development. His MATLAB code provided the basis
    for several transfer entropy methods and binning-related code.
- [George Datseris](https://github.com/Datseris) also ported KSG1 and KSG2 mutual
    information estimators to Neighborhood.jl.
- [Bjarte Hannisdal](https://github.com/bhannis) provided tutorials for mutual information.
- Tor Einar Møller contributed to cross-mapping methods in initial stages of development.

Many individuals has contributed code to other packages
in the [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/) ecosystem which
we use here. Contributors are listed in the respective GitHub repos and webpages.
