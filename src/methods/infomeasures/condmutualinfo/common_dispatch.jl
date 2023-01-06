using ComplexityMeasures: encode, decode
using ComplexityMeasures: RectangularBinEncoding
using ComplexityMeasures: FixedRectangularBinning, ValueHistogram
using StateSpaceSets: Dataset, AbstractDataset
using StaticArrays: @MVector

###########################################################################################
# Generic dispatch.
# ----------------
# Many `ProbabilitiesEstimator`s and `DifferentialEntropyEstimator`s can be directly
# applied to compute mutual information.
###########################################################################################
function estimate(measure::CMI{E, <:CMIMI2}, est::MutualInformationEstimator, x, y, z) where {E}
    X = Dataset(x)
    Y = Dataset(y)
    Z = Dataset(z)
    YZ = Dataset(Y, Z)
    m = measure.definition.measure
    return mutualinfo(m, est, X, YZ) + mutualinfo(m, est, X, Y)
end

function estimate(measure::CMI{E, D}, est, args...) where {E, D}
    estimate(default_measure(measure, est), est, args...)
    # t = typeof(est)
    # @warn """MutualInformationEstimator with definition $D can't be used with mutual
    #     information estimators such as $t.
    #     Please provide a `MutualInformationDefinition` that is compatible with $t and $E.
    #     See online documentation for an overview.\
    #     """
    #
end


###########################################################################################
# Specialized dispatch.
# ---------------------
# But not all `ProbabilitiesEstimator`s can be directly applied, because there is no
# straight-forward extension to CMI. Circumventing this is easy, however.
# We can arbitrarily define custom outcome spaces based on any
# ProbabilitiesEstimator`, and compute CMI based on probabilities estimated over
# these outcome spaces.
# Below, we specify which methods get special treatment. The key ingredient here
# is `marginal_probabilities`, which defines how marginal probabilities are computed
# for a given measure with a certain definition.
###########################################################################################

const WellDefinedCMIH4ProbEsts{m, D} = Union{
    SymbolicPermutation{m},
    ValueHistogram{<:FixedRectangularBinning{D}},
    Dispersion
} where {m, D}

function estimate(measure::CMI, est::WellDefinedCMIH4ProbEsts{m, D},
        x, y, z) where {m, D}
    e = measure.e
    pXZ, pYZ, pXYZ, pZ = marginal_probabilities(measure, est, x, y, z)
    HXZ = entropy(e, pXZ)
    HYZ = entropy(e, pYZ)
    HXYZ = entropy(e, pXYZ)
    HZ = entropy(e, pZ)
    return HXZ + HYZ - HXYZ - HZ
end
