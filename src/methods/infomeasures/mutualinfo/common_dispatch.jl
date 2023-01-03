using ComplexityMeasures: RectangularBinEncoding, RectangularBinning, AbstractBinning
using ComplexityMeasures: PermutationProbabilitiesEstimator
using ComplexityMeasures: entropy, encode, decode, symbolize_for_dispersion
using StateSpaceSets: AbstractDataset, Dataset
using StateSpaceSets: dimension
using DelayEmbeddings: embed, genembed

export marginal_probabilities


###########################################################################################
# Generic dispatch.
# ----------------
# Many `ProbabilitiesEstimator`s and `DifferentialEntropyEstimator`s can be directly
# applied to compute mutual information.
###########################################################################################
function estimate(measure::MutualInformation{<:EntropyDefinition, <:MIH3},
        est::ProbOrDiffEst,
        x, y)
    e = measure.e
    X = Dataset(x); Y = Dataset(y); XY = Dataset(X, Y)
    hX = entropy(e, est, X)
    hY = entropy(e, est, Y)
    hXY = entropy(e, est, XY)
    return hX + hY - hXY
end

###########################################################################################
# Specialized dispatch.
# ---------------------
# Not all `ProbabilitiesEstimator`s can be directly applied, because there is no
# straight-forward extension to mutual information. Circumventing this is easy, however.
# We can arbitrarily define custom outcome spaces based on any
# `ProbabilitiesEstimator`, and compute mutual information based on probabilities
# estimated over these outcome spaces.
# Below, we specify which methods get special treatment. The key ingredient here
# is `marginal_probabilities`, which defines how marginal probabilities are computed
# for a given measure with a certain definition.
###########################################################################################

const WellDefinedMIH3Probs{m, D} = Union{
    SymbolicPermutation{m},
    ValueHistogram{<:FixedRectangularBinning{D}},
    Dispersion
} where {m, D}

# Generic for all H3-type definitions.
function estimate(measure::MutualInformation{<:EntropyDefinition, <:MIH3},
        est::WellDefinedMIH3Probs{m, D},
        x, y) where {m, D}

    pX, pY, pXY = marginal_probabilities(measure, est, x, y)
    e = measure.e
    HX = entropy(e, pX)
    HY = entropy(e, pY)
    HXY = entropy(e, pXY)
    return HX + HY - HXY
end

function estimate(measure::MITsallis{<:EntropyDefinition, MIDefinitionTsallisH3Martin},
        est::WellDefinedMIH3Probs{m, D},
        x, y) where {m, D}

    pX, pY, pXY = marginal_probabilities(measure, est, x, y)
    e = measure.e
    HX = entropy(e, pX)
    HY = entropy(e, pY)
    HXY = entropy(e, pXY)
    q = measure.e.q
    return HX + HY - (1 - q) * HX * HY - HXY
end
