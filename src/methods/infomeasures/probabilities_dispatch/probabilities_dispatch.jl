export marginal_probabilities

"""
    marginal_probabilities(measure::InformationMeasure, est, x₁, x₂, …) → px₁, px₂, …

Compute marginal probabilities required for the estimator `est` to compute the
given `measure` (e.g. [`CMIShannon`](@ref)), given input data `x₁`, `x₂`, etc.

All methods of this function must guarantee that the returned probabilities `px₁, px₂, …`
all have the same number of elements (i.e. the joint outcome space is well-defined, and
the marginal probabilities are simply obtained by marginalizing the joint probabilities).
"""
function marginal_probabilities end

############################################################################################
# The following files contain relevant dispatch for all discrete probability-based
# estimators. The files are (arbitrarily) grouped by estimator.
############################################################################################

include("value_histogram.jl")
include("symbolic_permutation.jl")
include("dispersion.jl")
