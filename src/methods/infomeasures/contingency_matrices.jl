export ContingencyMatrix
export marginal_probs

import Base: size, getindex, setindex
import ComplexityMeasures: probabilities, outcomes

# Contigency matrices (i.e normalized multidimensional histograms) for 2D and 3D.
# This could in principle be done with fewer lines of code using a generated function,
# but we only need a maximum of 3D matrices, so we just do it manually for now.
# This is *not* optimized for speed, but it *works*

"""
    ContingencyMatrix{T, N} <: AbstractArray{T, N}
    ContingencyMatrix(x, y, [z], ...) → c::ContingencyMatrix
    ContingencyMatrix(est::ProbabilitiesEstimator, x, y, [z], ...) → c::ContingencyMatrix

Estimate a contingency matrix (an `N`-dimensional empirical probability mass function,
taking naive relative frequency counts as probability estimates).

- If the input data `x` and `y` (and additional variables like `z`, if provided) are
    already discretized, a contingency matrix is computed directly from the data and
    the observed outcomes are taken as the unique values occurring in the input data.
    Outcomes are ordered according to the returned outcomes from `sort(unique(X))`, etc.

- If a [`ProbabilitiesEstimator`](@ref) is applied, then the input variables `x`, `y`, `z`,
    etc, are discretized separately according to `est`. The contingency matrix is then
    computed based on the observed encodings/symbols/outcomes for each variable.
    This variant enforces *the same* outcome space for all variables, and calls
    [`marginal_encodings`](@ref) to get the discretized values.

`ContingencyMatrix` is just a simple wrapper around `AbstractArray{T, N}` that
represents multidimensional normalized frequency counts, i.e. a multivariate
version of [`Probabilities`](@ref). If indexed with two indices `i, j`, then
it returns the joint probability for the `(i, j)`th outcome.

The function `probabilities(x::ContingencyMatrix, i::Int)` returns the marginal
probabilities for the `i`-th dimension.
The function `outcomes(x::ContingencyMatrix, i::Int)` returns the marginal
outcomes for the `i`-th dimension.

## Usage

Contingency matrices is used in the computation of discrete versions of the following
quantities:

- [`entropy_joint`](@ref).
- [`mutualinfo`](@ref).
- [`condmutualinfo`](@ref).
"""
struct ContingencyMatrix{T, N, O} <: AbstractArray{T, N}
    # We only need to keep track of the joint probabilities. Marginals are obtained through,
    # unsurprisingly, marginalization of the joint probabilities.
    joint::AbstractArray{T, N}
    marginals::NTuple{N, Probabilities}
    marginal_outcomes::NTuple{N, O}
end
Base.size(c::ContingencyMatrix) = size(c.joint)
Base.getindex(c::ContingencyMatrix, i) = getindex(c.joint, i)
Base.getindex(c::ContingencyMatrix, i::Int...) = getindex(c.joint, i...)
Base.setindex!(c::ContingencyMatrix, v, i) = setindex!(c.joint, v, i)

probabilities(c::ContingencyMatrix, i::Int) = c.marginals[i]
outcomes(c::ContingencyMatrix, i::Int) = c.marginal_outcomes[i]

# TODO: actually dispatch on joint frequency method for all methods below.
function ContingencyMatrix(X, Y)
    pX, ΩX = probabilities_and_outcomes(CountOccurrences(), X); lX = length(pX)
    pY, ΩY = probabilities_and_outcomes(CountOccurrences(), Y); lY = length(pY)
    p̂X = reshape(pX, lX, 1)
    p̂Y = reshape(pY, 1, lY)
    pXY = p̂X .* p̂Y
    return ContingencyMatrix(pXY, (pX, pY), (ΩX, ΩY))
end

function ContingencyMatrix(X, Y, Z)
    pX, ΩX = probabilities_and_outcomes(CountOccurrences(), X); lX = length(pX)
    pY, ΩY = probabilities_and_outcomes(CountOccurrences(), Y); lY = length(pY)
    pZ, ΩZ = probabilities_and_outcomes(CountOccurrences(), Z); lZ = length(pZ)
    p̂X = reshape(pX, lX, 1, 1)
    p̂Y = reshape(pY, 1, lY, 1)
    p̂Z = reshape(pZ, 1, 1, lZ)
    pXYZ = p̂X .* p̂Y .* p̂Z
    return ContingencyMatrix(pXYZ, (pX, pY, pZ), (ΩX, ΩY, ΩZ))
end

function ContingencyMatrix(est::ProbabilitiesEstimator, x...)
    # Enforce same outcome space for all variables.
    encodings = marginal_encodings(est, x...)
    return ContingencyMatrix(encodings...)
end


# The following commented-out code below is equivalent to theabove, but muuuch faster.
# I keep the comments here for, so when I revisit this, I understand *why* it works.
# This will be moved into some sort of tutorial or doc example at some point.
# function ContingencyMatrix(X, Y)
    # Ωx = sort(unique(X))
    # Ωy = sort(unique(Y))
    # N = length(X) * length(Y)
    # pXY = zeros(length(Ωx), length(Ωy)) # enumerates the possible states
    # for (i, ωyᵢ) in enumerate(Ωy)
    #     for (j, ωxᵢ) in enumerate(Ωx)
    #         pXY[j, i] = count_occurrences(ωxᵢ, ωyᵢ, X, Y)
    #     end
    # end
    # return ContingencyMatrix(pXY / N)
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

# function count_occurrences(ωxᵢ, ωyᵢ, X, Y)
#     n = 0.0
#     for x in X
#         for y in Y
#             if x == ωxᵢ && y == ωyᵢ
#                 n += 1
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
