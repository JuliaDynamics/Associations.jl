# [Syntax overview](@id syntax_overview)

This is an overview over the low-level functionality in the CausalityTools.jl package and its subpackages.

## Discretization

- [`RectangularBinning`](@ref). Instructions for using rectangular partitions.
- [`TriangulationBinning`](@ref). Instructions for using triangulated partitions.

## Transfer operator estimation

- [`transferoperator(points, binning_scheme::RectangularBinning)`](@ref)

- [`transferoperator(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)`](@ref)

- [`transferoperator(pts, ϵ::TriangulationBinning, simplex_intersection_type::ApproximateIntersection)`](@ref)

## Invariant measure estimation

- [`invariantmeasure(points, binning_scheme::RectangularBinning)`](@ref)

- [`invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)`](@ref)

- [`invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ApproximateIntersection)`](@ref)

## Transfer entropy estimation

There are two estimators that compute transfer entropy by rectangular partitions. 

- [`transferentropy(pts, vars::TEVars, ϵ::RectangularBinning, estimator::VisitationFrequency; b = 2)`](@ref)

- [`transferentropy(pts, vars::TEVars, ϵ::RectangularBinning, estimator::TransferOperatorGrid; b = 2)`](@ref)

To compute transfer entropy over triangulated partitions, the invariant measure over the 
triangulation must be precomputed, using either 

- [`invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ExactIntersection)`](@ref) or 
- [`invariantmeasure(pts, ϵ::TriangulationBinning, simplex_intersection_type::ApproximateIntersection)`](@ref)).

Then we can superimpose rectangular grid over the triangulation and compute the transfer entropy.

- [`transferentropy(μ::AbstractTriangulationInvariantMeasure, vars::TEVars, binning_scheme::RectangularBinning; estimator = VisitationFrequency(), n::Int = 10000, b = 2)`](@ref).

## Cross mappings

- [`crossmap`](@ref).
- [`convergentcrossmap`](@ref).

## Joint distance distribution

- [`joint_distance_distribution`](@ref)
