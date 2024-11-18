using Neighborhood: Euclidean, KDTree, NeighborNumber, Theiler
using Neighborhood: bulksearch, inrangecount

export AzadkiaChatterjeeCoefficient

"""
    AzadkiaChatterjeeCoefficient <: AssociationMeasure
    AzadkiaChatterjeeCoefficient(; theiler::Int = 0)

The Azadkia-Chatterjee coefficient [Azadkia2021](@cite) is a coefficient for pairwise and conditional 
association inspired by the Chatterjee-Dette-Siburg-Stoimenov coefficient [Chatterjee2021, Dette2013](@cite)
(see [`ChatterjeeCorrelation`](@ref)).

# Usage

- Use with [`association`](@ref) to compute the raw Azadkia-Chatterjee coefficient.
- Use with [`SurrogateAssociationTest`](@ref) to perform a surrogate test for significance 
    of a pairwise or conditional Azadkia-Chatterjee-type association ([example](@ref example_SurrogateAssociationTest_AzadkiaChatterjeeCoefficient)). 
    When using a surrogate test for significance, only the *first* input variable is shuffled 
    according to the given surrogate method.
- Use with [`LocalPermutationTest`](@ref) to perform a test of conditional independence 
    ([example](@ref example_LocalPermutationTest_AzadkiaChatterjeeCoefficient)).

## Description

The pairwise statistic is 

```math
T_n(Y, \\boldsymbol{Z} | \\boldsymbol{X}) = \\dfrac{\\sum_{i=1}^n \\left( \\min{(R_i, R_{M_{(i)}})} - \\min{(R_i, R_{N_{(i)}})} \\right) }{\\sum_{i=1}^n \\left(R_i - \\min{(R_i, R_{N_{(i)}})} \\right)}.
```

where ``R_i`` is the rank of the point ``Y_i`` among all ``Y_i``s, and ``M_{(i)}`` and ``N_{(i)}`` are indices 
of nearest neighbors of the points ``\\boldsymbol{X}_i`` and ``(\\boldsymbol{X}_i, \\boldsymbol{Z}_i)``, respectively
(given appropriately constructed marginal spaces). The `theiler` keyword argument is an integer controlling the 
number of nearest neighbors to exclude during neighbor searches. The Theiler window defaults to `0`, which excludes 
self-neighbors, and is the only option considered in [Azadkia2021](@citet).

In the case where ``\\boldsymbol{X}`` has no components (i.e. we're not conditioning), we also consider 
``L_i`` as the number of ``j`` such that ``Y_j \\geq Y_i``. The measure is then defined as 

```math
T_n(Y, \\boldsymbol{Z}) = 
\\dfrac{\\sum_{i=1}^n \\left( n \\min{(R_i, R_{M_{(i)}})} - L_i^2 \\right) }{\\sum_{i=1}^n \\left( L_i (n - L_i) \\right)}.
```

The value of the coefficient is on `[0, 1]` when the number of samples goes to `∞`, but is not restricted
to this interval in practice. 

## Input data

If the input data contain duplicate points, consider adding a small magnitude of noise to the input 
data. Otherwise, errors will occur when locating nearest neighbors.

# Estimation

- [Example 1](@ref example_AzadkiaChatterjeeCoefficient). Estimating the Azadkia-Chatterjee
    coefficient to quantify associations for a chain of unidirectionally coupled variables,
    showcasing both pairwise and conditional associations.
- [Example 2](@ref example_SurrogateAssociationTest_AzadkiaChatterjeeCoefficient). Using 
    [`SurrogateAssociationTest`](@ref) in combination with the Azadkia-Chatterjee coefficient
    to quantify significance of pairwise and conditional associations.
- [Example 3](@ref example_LocalPermutationTest_AzadkiaChatterjeeCoefficient). 
    Using [`LocalPermutationTest`](@ref) in combination with the Azadkia-Chatterjee coefficient
    to perform a test for conditional independence.
"""
struct AzadkiaChatterjeeCoefficient{V, L, R, MI, NI, RM, RN} <: AssociationMeasure
    sorted_y::V
    Lᵢs::L
    Rᵢs::R
    Mᵢs::MI
    Nᵢs::NI
    RMᵢs::RM
    RNᵢs::RN
    theiler::Int
end

function Base.show(io::IO, m::AzadkiaChatterjeeCoefficient)
    theiler = m.theiler
    print(io, "AzadkiaChatterjeeCoefficient(; theiler = $theiler)")
end


min_inputs_vars(::AzadkiaChatterjeeCoefficient) = 2
max_inputs_vars(::AzadkiaChatterjeeCoefficient) = 3

function AzadkiaChatterjeeCoefficient(; theiler::I = 0) where I <: Integer
    return AzadkiaChatterjeeCoefficient{Nothing, Nothing, Nothing, Nothing, Nothing, Nothing, Nothing}(
        nothing, nothing, nothing, nothing, nothing, nothing, nothing, theiler,
    )
end

function AzadkiaChatterjeeCoefficient(inputs...; theiler::Integer = 0)
    @assert 2 ≤ length(inputs) ≤ 3
    ns = length.(inputs)
    @assert allequal(ns)
    m = length(inputs)
    n = length(first(inputs))
    sorted_ys = sort(first(inputs))

    r = zeros(Int, n)

    if m == 2
        r = zeros(Int, n)
        Mᵢs = zeros(Int, n)
        Nᵢs = nothing
        RMᵢs = zeros(Int, n)
        RNᵢs = nothing
        𝓁 =  zeros(Int, n)
    else # m == 3
        r = zeros(Int, n)
        Mᵢs = zeros(Int, n)
        Nᵢs = zeros(Int, n)
        RMᵢs = zeros(Int, n)
        RNᵢs = zeros(Int, n)
        𝓁 = nothing
    end

    return AzadkiaChatterjeeCoefficient(sorted_ys, 𝓁, r, Mᵢs, Nᵢs, RMᵢs, RNᵢs, theiler)
end

# For non-preallocated version, simply create a pre-allocated copy of the definition
# and call the code that uses the pre-allocated definition.
function association(m::AzadkiaChatterjeeCoefficient{Nothing}, inputs...)
    m_preallocated = AzadkiaChatterjeeCoefficient(inputs...; theiler = m.theiler)
    return association(m_preallocated, inputs...)
end

# todo: add note about adding some slight noise of there are identical data points.
# todo: separate method for two variables
# note: the input order follows the notation in AZADKIA & chatterjee.
# The value of Tₙ is in the limit of n → ∞ guaranteed to be on [0, 1], but will not be so in practice.
function association(m::AzadkiaChatterjeeCoefficient{<:AbstractVector}, 
        y::Union{AbstractVector{<:Real}, AbstractStateSpaceSet{1}}, 
        z::VectorOrStateSpaceSet, 
        x::VectorOrStateSpaceSet) 
    @assert length(y) == length(z) == length(x)
    n = length(x); @assert n ≥ 2

    metric = Euclidean()

    X = StateSpaceSet(x)
    XZ = StateSpaceSet(x, z)
    tree_X = KDTree(X, metric)
    tree_XZ = KDTree(XZ, metric)

    W = Theiler(m.theiler)

    Nᵢs = first.(bulkisearch(tree_X, X, NeighborNumber(1), W))
    Mᵢs = first.(bulkisearch(tree_XZ, XZ, NeighborNumber(1), W))
    loopcount_rs!(m.Rᵢs, m.sorted_y, y) # from chatterjee.jl

    # presort
    m.RMᵢs .= m.Rᵢs[Mᵢs]
    m.RNᵢs .= m.Rᵢs[Nᵢs]

    num = 0.0
    den = 0.0
    for i = 1:n
        Rᵢ = m.Rᵢs[i]
        RMᵢ = m.RMᵢs[i]
        RNᵢ = m.RNᵢs[i]
        mn = min(Rᵢ, RNᵢ)
        num += (min(Rᵢ, RMᵢ) - mn)
        den += Rᵢ - mn
    end
    Tₙ = num / den
    return Tₙ
end

function association(m::AzadkiaChatterjeeCoefficient{<:AbstractVector}, y::AbstractVector{<:Real}, z::VectorOrStateSpaceSet) 
    @assert length(y) == length(z)
    n = length(y); @assert n ≥ 2

    metric = Euclidean()

    Z = StateSpaceSet(z)
    tree_Z = KDTree(Z, metric)
    W = Theiler(m.theiler)
    Mᵢs = first.(bulkisearch(tree_Z, Z, NeighborNumber(1), W))

    loopcount_rs_ℓs!(m.Rᵢs, m.Lᵢs, m.sorted_y, y) # from chatterjee.jl
    
    # presort
    m.RMᵢs .= m.Rᵢs[Mᵢs]

    num = 0.0
    den = 0.0
    for i = 1:n
        Lᵢ = m.Lᵢs[i]
        Rᵢ = m.Rᵢs[i]
        RMᵢ = m.RMᵢs[i]
        num += n * min(Rᵢ, RMᵢ) - Lᵢ^2
        den += Lᵢ * (n - Lᵢ)
    end
    Tₙ = num / den
    return Tₙ
end