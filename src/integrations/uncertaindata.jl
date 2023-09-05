# The uncertainty handling framework in this file will be added
# as part of a 2.X release. Can be ignored for now.

import UncertainData:
    resample,
    UncertainStateSpaceSet,
    UncertainIndexStateSpaceSet,
    UncertainValueStateSpaceSet,
    UncertainIndexValueStateSpaceSet
import .s_measure
import .jdd
import .transferentropy; export transferentropy
import .crossmap
import .ccm
import .predictive_asymmetry
import HypothesisTests.OneSampleTTest

##################################################
# Basic resampling for `UncertainStateSpaceSet`s
##################################################
const UT = Union{UncertainValueStateSpaceSet, UncertainIndexStateSpaceSet, UncertainStateSpaceSet}

s_measure(s::UT, t::UT, args...; kwargs...) =
    s_measure(resample(s), resample(t), args...; kwargs...)

jdd(s::UT, t::UT; kwargs...) =
    jdd(resample(s), resample(t); kwargs...)

jdd(test::OneSampleTTest, s::UT, t::UT; kwargs...) =
    jdd(test, resample(s), resample(t); kwargs...)

mutualinfo(s::UT, t::UT, method; kwargs...) =
    mutualinfo(resample(s), resample(t), method; kwargs...)

info_methods = [
    :VisitationFrequency, :TransferOperator,
    :OrdinalPatterns, :AmplitudeAwareOrdinalPatterns, :WeightedOrdinalPatterns,
    :NaiveKernel,
    :Kraskov,
    :Kraskov1,
    :Kraskov2,
    :KozachenkoLeonenko,
    :Hilbert,
    :TimeScaleMODWT
]

for method in info_methods
    @eval transferentropy($(method), s::UT, t::UT; kwargs...) =
        transferentropy(method, resample(s), resample(t); kwargs...)

    @eval transferentropy($(method), s::UT, t::UT, c::UT; kwargs...) =
        transferentropy(method, resample(s), resample(t), resample(c); kwargs...)
end

# transferentropy(s::UT, t::UT, method; kwargs...) =
#     transferentropy(resample(s), resample(t), method; kwargs...)

# transferentropy(s::UT, t::UT, c::UT, method; kwargs...) =
#     transferentropy(resample(s), resample(t), resample(c), method; kwargs...)

predictive_asymmetry(method, s::UT, t::UT; kwargs...) =
    predictive_asymmetry(method, resample(s), resample(t); kwargs...)

predictive_asymmetry(method, s::UT, t::UT, c::UT; kwargs...) =
    predictive_asymmetry(method, resample(s), resample(t), resample(c); kwargs...)

crossmap(s::UT, t::UT, args...; kwargs...) =
    crossmap(resample(s), resample(t), args...; kwargs...)

ccm(s::UT, t::UT, args...; kwargs...) =
    ccm(resample(s), resample(t), args...; kwargs...)

##########################################################################
# Basic resampling for `UncertainIndexValueStateSpaceSet` (no constraints)
##########################################################################
const UIVD = UncertainIndexValueStateSpaceSet

# TODO: warn about potential index reversals?
#
# function warn_about_sampling(s::V, t::W)
#     if s isa UIVD
#         @warn "`s` isa UncertainIndexValueStateSpaceSet. Index reversals may occur. Consider constrained resampling."
#     end

#     if t isa UIVD
#         @warn "`t` isa UncertainIndexValueStateSpaceSet. Index reversals may occur. Consider constrained resampling."
#     end
# end

# function warn_about_sampling(s::V, t::W, c::X)
#     warn_about_sampling(s, t)
#     if c isa UIVD
#         @warn "`c` isa UncertainIndexValueStateSpaceSet. Index reversals may occur. Consider constrained resampling."
#     end
# end

s_measure(s::UIVD, t::UIVD, args...; kwargs...) =
    s_measure(resample(s.values), resample(t.values), args...; kwargs...)

jdd(s::UIVD, t::UIVD; kwargs...) =
    jdd(resample(s.values), resample(t.values); kwargs...)

jdd(test::OneSampleTTest, s::UIVD, t::UIVD; kwargs...) =
    jdd(test, resample(s), resample(t); kwargs...)

mutualinfo(method, s::UIVD, t::UIVD; kwargs...) =
    mutualinfo(method, resample(s.values), resample(t.values); kwargs...)

transferentropy(method, s::UIVD, t::UIVD; kwargs...) =
    transferentropy(method, resample(s.values), resample(t.values); kwargs...)

transferentropy(method, s::UIVD, t::UIVD, c::UIVD; kwargs...) =
    transferentropy(method, resample(s.values), resample(t.values), resample(c.values); kwargs...)

predictive_asymmetry(method, s::UIVD, t::UIVD; kwargs...) =
    predictive_asymmetry(method, resample(s.values), resample(t.values); kwargs...)

predictive_asymmetry(s::UIVD, t::UIVD, c::UIVD, method; kwargs...) =
    predictive_asymmetry(resample(s.values), resample(t.values), resample(c.values), method; kwargs...)

crossmap(s::UIVD, t::UIVD, args...; kwargs...) =
    crossmap(resample(s.values), resample(t.values), args...; kwargs...)

ccm(s::UIVD, t::UIVD, args...; kwargs...) =
    ccm(resample(s.values), resample(t.values), args...; kwargs...)
