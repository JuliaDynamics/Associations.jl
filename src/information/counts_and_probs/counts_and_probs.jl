export Discretization

# The type parameter `N` indicates the number of input datasets to be discretized.
"""
    Discretization

The supertype of all discretization schemes.

## Description 

A fundamental operation when computing multivariate information measures from data 
is *discretization*. There are many ways of discretizing multiple input datasets. We 
offer two main ways of doing so.

## Concrete implementations

- [`CodifyVariables`](@ref)
- [`CodifyPoints`](@ref)
"""
abstract type Discretization{N} end


# Concrete ways of encoding multivariate data, each defined as a type.
include("encoding/codify_points.jl")
include("encoding/codify_variables.jl")

# Counting and probabilities (contingency tables and probabilities for multivariate data)
include("counts/counts.jl")
include("probabilities/probabilities.jl")

# Estimating counts/probabilities from multiple input datasets.
include("estimation_per_point.jl")
include("estimation_per_variable.jl")
