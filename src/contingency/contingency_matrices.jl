import Base: size, getindex, setindex
import ComplexityMeasures: probabilities, outcomes
using StatsBase: levelsmap

export ContingencyMatrix
export probabilities, outcomes, frequencies
export contingency_matrix

# This should be done so that the following can be added to the docs, but for now,
# we only need the table, not information about variables.
# Let `c` be a 2-dimensional `ContingencyMatrix` constructed from input data `x` and `y`.
# Let `Ωx = unique(x)` and `Ωy = unique(y)`. `c[i, j]` then corresponds to the
# outcome `(unique(Ωx)[i], unique(Ωy)[j]`). The generalization to higher dimensions is
# straight-forward.

# This is not optimized for speed, but it *works*
"""
    ContingencyMatrix{T, N} <: Probabilities{T, N}
    ContingencyMatrix(frequencies::AbstractArray{Int, N})

A contingency matrix is essentially a multivariate analogue of [`Probabilities`](@ref)
that also keep track of raw frequencies.

The contingency matrix can be constructed directyly from an `N`-dimensional `frequencies`
array. Alternatively, the [`contingency_matrix`](@ref) function performs counting for
you; this works on both raw categorical data, or by first discretizing data using a
a [`ProbabilitiesEstimator`](@ref).

## Description

A `ContingencyMatrix` `c` is just a simple wrapper around around `AbstractArray{T, N}`.
Indexing `c` with multiple indices `i, j, …` returns the `(i, j, …)`th
element of the empirical probability mass function (pmf). The following convencience
methods are defined:

- `frequencies(c; dims)` returns the multivariate raw counts along the given `dims
    (default to all available dimensions).
- `probabilities(c; dims)` returns a multidimensional empirical
    probability mass function (pmf) along the given `dims` (defaults to all available
    dimensions), i.e. the normalized counts.
- `probabilities(c, i::Int)` returns the marginal probabilities for the `i`-th dimension.
- `outcomes(c, i::Int)` returns the marginal outcomes for the `i`-th dimension.

# Ordering

The ordering of outcomes are internally consistent, but we make no promise on the
ordering of outcomes relative to the input data. This means that if your input
data are `x = rand(["yes", "no"], 100); y = rand(["small", "medium", "large"], 100)`,
you'll get a 2-by-3 contingency matrix, but there currently no easy way to
determine which outcome the i-j-th row/column of this matrix corresponds to.

Since [`ContingencyMatrix`](@ref) is intended for use in information theoretic methods
that don't care about ordering, as long as the ordering is internally consistent,
this is not an issue for practical applications in this package.
This may change in future releases.

## Usage

Contingency matrices is used in the computation of discrete versions of the following
quantities:

- [`entropy_joint`](@ref).
- [`mutualinfo`](@ref).
- [`condmutualinfo`](@ref).
"""
struct ContingencyMatrix{T, N, I} <: AbstractArray{T, N}
    # We only need to keep track of the joint probabilities. Marginals are obtained through,
    # unsurprisingly, marginalization of the joint probabilities.
    probs::AbstractArray{T, N}
    freqs::AbstractArray{I, N}
end

function ContingencyMatrix(freqs::AbstractArray{Int, N}) where N
    probs = freqs ./ sum(freqs)
    return ContingencyMatrix(probs, freqs)
end

Base.size(c::ContingencyMatrix) = size(c.probs)
Base.getindex(c::ContingencyMatrix, i) = getindex(c.probs, i)
Base.getindex(c::ContingencyMatrix, i::Int...) = getindex(c.probs, i...)
Base.setindex!(c::ContingencyMatrix, v, i) = setindex!(c.probs, v, i)

function frequencies(c::ContingencyMatrix; dims = 1:ndims(c))
    alldims = 1:ndims(c)
    reduce_dims = setdiff(alldims, dims)
    marginal = dropdims(sum(c.freqs, dims = reduce_dims), dims = (reduce_dims...,))
    return marginal
end

function probabilities(c::ContingencyMatrix; dims = 1:ndims(c))
    alldims = 1:ndims(c)
    reduce_dims = (setdiff(alldims, dims)...,)
    marginal = dropdims(sum(c.probs, dims = reduce_dims), dims = reduce_dims)
    return Probabilities(marginal)
end

"""
    contingency_matrix(x, y, [z, ...]) → c::ContingencyMatrix
    contingency_matrix(est::ProbabilitiesEstimator, x, y, [z, ...]) → c::ContingencyMatrix

Estimate a multidimensional contingency matrix `c` from input data `x, y, …`, where the
input data can be of any and different types, as long as `length(x) == length(y) == …`.

For already discretized data, use the first method. For continuous data, you want to
discretize the data before computing the contingency table. You can do
this manually and then use the first method. Alternatively, you can provide a
[`ProbabilitiesEstimator`](@ref) as the first
argument to the constructor. Then the input variables `x, y, …` are discretized
*separately* according to `est` (*enforcing the same outcome space for all variables*),
by calling [`marginal_encodings`](@ref).
"""
function contingency_matrix end

function contingency_matrix(est::ProbabilitiesEstimator, x...)
    # Enforce identical outcome space for all variables.
    encodings = marginal_encodings(est, x...)
    return contingency_matrix(encodings...)
end

# For this to work generically, we must map unique elements to integers.
function contingency_matrix(x...)
    Ls = length.(x);
    if !all(Ls .== maximum(Ls))
        throw(ArgumentError("Input data must have equal lengths. Got lengths $Ls."))
    end

    # The frequency matrix dimensions are dictated by the number of unique occurring values
    Ωs = unique.(x)
    matrix_dims = length.(Ωs);

    # Get marginal probabilities and outcomes
    #pΩs = [probabilities_and_outcomes(CountOccurrences(), xᵢ) for xᵢ in x]
    freqs, lmaps = freqtable_equallength(matrix_dims, x...)

    # TODO: Inverse map from integer-encoded outcomes to the original outcomes.
    # marginal_outcomes = [map(k -> lmap[k], last(pΩ)) for (pΩ, lmap) in zip(pΩs, lmaps)]
    probs = freqs ./ maximum(Ls)
    return ContingencyMatrix(
        probs,
        freqs,
    )
end


function freqtable_equallength(matrix_dims, x...)
    # Map the input data to integers. This ensures compatibility with *any* input type.
    # Then, we can simply create a joint `Dataset{length(x), Int}` and use its elements
    # as `CartesianIndex`es to update counts.
    lvl = tolevels.(x)
    levels = (first(l) for l in lvl)
    lmaps = [last(l) for l in lvl]
    X = Dataset(levels...)

    freqs = zeros(Int, matrix_dims)
    for ix in to_cartesian(sort(X.data)) # sorted matrix access should be faster.
        freqs[ix] += 1
    end
    return freqs, lmaps
end

function to_cartesian(x)
    (CartesianIndex.(xᵢ...) for xᵢ in x)
end

"""
    tolevels!(levels, x) → levels, dict
    tolevels(x) → levels, dict

Apply the bijective map ``f : \\mathcal{Q} \\to \\mathbb{N}^+`` to each `x[i]` and store
the result in `levels[i]`, where `levels` is a pre-allocated integer vector such that
`length(x) == length(levels)`.

``\\mathcal{Q}`` can be any space, and each ``q \\in \\mathcal{Q}`` is mapped to a unique
integer  in the range `1, 2, …, length(unique(x))`. This is useful for integer-encoding
categorical data such as strings, or other complex data structures.

The single-argument method allocated a `levels` vector internally.

`dict` gives the inverse mapping.
"""
function tolevels!(levels, x)
    @assert length(levels) == length(x)
    lmap = levelsmap(x)
    for i in eachindex(x)
        levels[i] = lmap[x[i]]
    end
    return levels, lmap
end

function tolevels(x)
    lmap = levelsmap(x)
    levels = zeros(Int, length(x))
    for i in eachindex(x)
        levels[i] = lmap[x[i]]
    end
    return levels, lmap
end

include("Contingency.jl")


# The following commented-out code below is equivalent to theabove, but muuuch faster.
# I keep the comments here for, so when I revisit this, I understand *why* it works.
# This will be moved into some sort of tutorial or doc example at some point.


# TODO: actually dispatch on joint frequency method for all methods below.
# function ContingencyMatrix(X, Y)
#     pX, ΩX = probabilities_and_outcomes(CountOccurrences(), X); lX = length(pX)
#     pY, ΩY = probabilities_and_outcomes(CountOccurrences(), Y); lY = length(pY)
#     p̂X = reshape(pX, lX, 1)
#     p̂Y = reshape(pY, 1, lY)
#     pXY = p̂X .* p̂Y
#     return ContingencyMatrix(pXY, pXY, (pX, pY), (ΩX, ΩY))
# end

# function ContingencyMatrix(X, Y, Z)
#     pX, ΩX = probabilities_and_outcomes(CountOccurrences(), X); lX = length(pX)
#     pY, ΩY = probabilities_and_outcomes(CountOccurrences(), Y); lY = length(pY)
#     pZ, ΩZ = probabilities_and_outcomes(CountOccurrences(), Z); lZ = length(pZ)
#     p̂X = reshape(pX, lX, 1, 1)
#     p̂Y = reshape(pY, 1, lY, 1)
#     p̂Z = reshape(pZ, 1, 1, lZ)
#     pXYZ = p̂X .* p̂Y .* p̂Z
#     return ContingencyMatrix(pXYZ, (pX, pY, pZ), (ΩX, ΩY, ΩZ))
# end






#function ContingencyMatrix(X, Y, Z)
    # Ωx = sort(unique(X))
    # Ωy = sort(unique(Y))
    # Ωz = sort(unique(Z))
    # N = length(X)*length(Y)*length(Z)
    # pXYZ = zeros(length(Ωx), length(Ωy), length(Ωz)) # enumerates the states
    # for (h, ωzᵢ) in enumerate(Ωz)
    #     for (j, ωyᵢ) in enumerate(Ωy)
    #         for (i, ωxᵢ) in enumerate(Ωx)
    #             pXYZ[i, j, h] = count_occurrences(ωxᵢ, ωyᵢ, ωzᵢ, X, Y, Z)
    #         end
    #     end
    # end
    #return ContingencyMatrix(pXYZ / N)

    # The following is equivalent to the commented-out code above, but muuuch faster.
     # I keep the above code, so when I revisit this, I understand *why* it works.
#     pX = probabilities(CountOccurrences(), X); lX = length(pX)
#     pY = probabilities(CountOccurrences(), Y); lY = length(pY)
#     pZ = probabilities(CountOccurrences(), Z); lZ = length(pZ)

#     # # Reshape explicitly for 3D case to work.
#     p̂X = reshape(pX, lX, 1, 1)
#     p̂Y = reshape(pY, 1, lY, 1)
#     p̂Z = reshape(pZ, 1, 1, lZ)
#     pXYZ = p̂X .* p̂Y .* p̂Z

#     return pXYZ
# end

# function count_occurrences(ωxᵢ, ωyᵢ, ωzᵢ, X, Y, Z)
#     n = 0.0
#     for x in X
#         for y in Y
#             for z in Z
#                 if x == ωxᵢ && y == ωyᵢ && z == ωzᵢ
#                     n += 1
#                 end
#             end
#         end
#     end
#     return n
# end

# function count_occurrences(total_count, ωxᵢ, ωyᵢ, X, Y)
#     n = 0.0
#     for x in X
#         for y in Y
#             if x == ωxᵢ && y == ωyᵢ
#                 n += 1
#                 total_count[] += 1
#             end
#         end
#     end
#     return n
# end


###########################################################################################
# The commented-out code below will be part of a later release. It is based on the
# fact that the way in which joint probabilities are estimated using the
# counts greatly affects the measure that is estimated based on the contingency table,
# like for mutual information (Fernandes et al., 2010)
###########################################################################################

# Exaplanation of the algorithm
# First, we compute the counts of `fᵢⱼ` for each outcome `(Ωxᵢ, Ωyᵢ)`
# by counting how many times `X == Ωxᵢ` and `y == Ωyᵢ` occur simultaneously.

# Then, the joint probabilities are computed according to `method`, which controls
# which prior assumptions about the joint relationship is assumed
# (Fernandes et al., 2010)[^Fernandes2010].

# [^Fernandes2010]: Fernandes, A. D., & Gloor, G. B. (2010). Mutual information is
#     critically dependent on prior assumptions: would the correct estimate of mutual
#     information please identify itself?. Bioinformatics, 26(9), 1135-1139.

# """
#     JointFrequencyRelationship

# The supertype for all ways of computing joint frequency relationships for
# a [`ContingencyMatrix`](@ref).

# ## Motivation

# As discussed in Fernandes et al. (2010), the mutual information estimated based on
# contingency matrices *strongly* depends on how joint probabilities are estimated
# for the contingency matrices.
# They present multiple alternatives, based of different prior assuptions of the data,
# which we here implement as subtypes:

# - [`NaiveJointFrequency`](@ref).
# - [`MultiplicativeJointFrequency`](@ref).

# [^Fernandes2010]: Fernandes, A. D., & Gloor, G. B. (2010). Mutual information is
#     critically dependent on prior assumptions: would the correct estimate of mutual
#     information please identify itself?. Bioinformatics, 26(9), 1135-1139.
# """
# abstract type JointFrequencyRelationship end

# """
#     NaiveJointFrequency <: JointFrequencyRelationship
#     NaiveJointFrequency()

# Takes the joint frequencies naively as the observed sample counts, divided by the
# total number of samples.

# [^Fernandes2010]: Fernandes, A. D., & Gloor, G. B. (2010). Mutual information is
#     critically dependent on prior assumptions: would the correct estimate of mutual
#     information please identify itself?. Bioinformatics, 26(9), 1135-1139.
# """
# struct NaiveJointFrequency <: JointFrequencyRelationship end

# """
#     MultiplicativeJointFrequency  <: JointFrequencyRelationship
#     MultiplicativeJointFrequency ()

# Represents the assumption of independence of contingency table frequencies,
# such that the joint probabilities are given by the outer product of
# the marginal frequencies (e.g. row-sums and column sums in the 2D case)[^Fernandes2010].

# [^Fernandes2010]: Fernandes, A. D., & Gloor, G. B. (2010). Mutual information is
#     critically dependent on prior assumptions: would the correct estimate of mutual
#     information please identify itself?. Bioinformatics, 26(9), 1135-1139.
# """
# struct MultiplicativeJointFrequency <: JointFrequencyRelationship end
