
# TODO: Name this something reasonable. EstimatinMethod is generic,
"""
    EstimationMethod

Some information measures, like mutual information or conditional mutual information,
may be computed in sevaral different ways. In CausalityTools, we call these
*compound* measures. For example, mutual information may be computed as a sum of
marginal entropy terms, or in terms of a KL-divergence.

`EstimationMethod` is the supertype of all compound estimation methods, which
in the case of information measures are methods of decomposing higher-level measures
into lower-level ones. Currently, subtypes are

- [`MI2`](@ref) (used to estimate conditional mutual information)
- [`H4`](@ref) (used to estimate conditional mutual information)
- [`H3`](@ref) (used to estimate mutual information)
"""
abstract type EstimationMethod end
struct MI2 <: EstimationMethod end
struct H3 <: EstimationMethod end
struct H4 <: EstimationMethod end

struct KLDiv <: EstimationMethod end

export MI2, H3, H4

# Things that will be eventually moved to Entropies.jl
include("various/probabilities.jl")
include("various/entropies.jl")

# The below belong in this package.
include("entropy_cross/entropy_cross.jl")
include("entropy_relative/entropy_relative.jl")
include("entropy_conditional/entropy_conditional.jl")
include("mutualinfo/mutualinfo.jl")
include("conditional_mutualinfo/conditional_mutualinfo.jl")
#include("transferentropy/transferentropy.jl")
#include("predictive_asymmetry/predictive_asymmetry.jl")

# TODO:
# Explain somewhere in the documentation that if `ProbabilitiesEstimators`
# are used, then the *discrete* version of whatever measure is computed.
# Otherwise, the *differential* entropy/mi/cmi/whatever is estimated.
# Some methods are *compound measures*, in the sense that they can be built
# from lower-level measures.
