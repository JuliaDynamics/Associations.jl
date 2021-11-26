
import UncertainData: 
    resample, 
    UncertainDataset, 
    UncertainIndexDataset, 
    UncertainValueDataset, 
    UncertainIndexValueDataset
import .s_measure
import .jdd
import .transferentropy; export transferentropy
import .crossmap
import .ccm
import .predictive_asymmetry
import HypothesisTests.OneSampleTTest

##################################################
# Basic resampling for `UncertainDataset`s
##################################################
const UT = Union{UncertainValueDataset, UncertainIndexDataset, UncertainDataset}

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
    :SymbolicPermutation, :SymbolicAmplitudeAwarePermutation, :SymbolicWeightedPermutation, 
    :NaiveKernel,
    :Kraskov,
    :Kraskov1,
    :Kraskov2,
    :KozachenkoLeonenko,
    :Hilbert,
    :TimeScaleMODWT
]

for method in info_methods 
    @eval transferentropy(s::UT, t::UT, $(method); kwargs...) =
        transferentropy(resample(s), resample(t), method; kwargs...)

    @eval transferentropy(s::UT, t::UT, c::UT, $(method); kwargs...) =
        transferentropy(resample(s), resample(t), resample(c), method; kwargs...)
end

# transferentropy(s::UT, t::UT, method; kwargs...) =
#     transferentropy(resample(s), resample(t), method; kwargs...)

# transferentropy(s::UT, t::UT, c::UT, method; kwargs...) =
#     transferentropy(resample(s), resample(t), resample(c), method; kwargs...)

predictive_asymmetry(s::UT, t::UT, method; kwargs...) =
    predictive_asymmetry(resample(s), resample(t), method; kwargs...)

predictive_asymmetry(s::UT, t::UT, c::UT, method; kwargs...) =
    predictive_asymmetry(resample(s), resample(t), resample(c), method; kwargs...)

crossmap(s::UT, t::UT, args...; kwargs...) = 
    crossmap(resample(s), resample(t), args...; kwargs...)

ccm(s::UT, t::UT, args...; kwargs...) = 
    ccm(resample(s), resample(t), args...; kwargs...)

##########################################################################
# Basic resampling for `UncertainIndexValueDataset` (no constraints)
##########################################################################
const UIVD = UncertainIndexValueDataset

# TODO: warn about potential index reversals?
#
# function warn_about_sampling(s::V, t::W) 
#     if s isa UIVD
#         @warn "`s` isa UncertainIndexValueDataset. Index reversals may occur. Consider constrained resampling."
#     end

#     if t isa UIVD
#         @warn "`t` isa UncertainIndexValueDataset. Index reversals may occur. Consider constrained resampling."
#     end
# end

# function warn_about_sampling(s::V, t::W, c::X)
#     warn_about_sampling(s, t)
#     if c isa UIVD
#         @warn "`c` isa UncertainIndexValueDataset. Index reversals may occur. Consider constrained resampling."
#     end
# end

s_measure(s::UIVD, t::UIVD, args...; kwargs...) =
    s_measure(resample(s.values), resample(t.values), args...; kwargs...)

jdd(s::UIVD, t::UIVD; kwargs...) =
    jdd(resample(s.values), resample(t.values); kwargs...)

jdd(test::OneSampleTTest, s::UIVD, t::UIVD; kwargs...) =
    jdd(test, resample(s), resample(t); kwargs...)
    
mutualinfo(s::UIVD, t::UIVD, method; kwargs...) =
    mutualinfo(resample(s.values), resample(t.values), method; kwargs...)

transferentropy(s::UIVD, t::UIVD, method; kwargs...) =
    transferentropy(resample(s.values), resample(t.values), method; kwargs...)

transferentropy(s::UIVD, t::UIVD, c::UIVD, method; kwargs...) =
    transferentropy(resample(s.values), resample(t.values), resample(c.values), method; kwargs...)

predictive_asymmetry(s::UIVD, t::UIVD, method; kwargs...) =
    predictive_asymmetry(resample(s.values), resample(t.values), method; kwargs...)

predictive_asymmetry(s::UIVD, t::UIVD, c::UIVD, method; kwargs...) =
    predictive_asymmetry(resample(s.values), resample(t.values), resample(c.values), method; kwargs...)

crossmap(s::UIVD, t::UIVD, args...; kwargs...) = 
    crossmap(resample(s.values), resample(t.values), args...; kwargs...)

ccm(s::UIVD, t::UIVD, args...; kwargs...) = 
    ccm(resample(s.values), resample(t.values), args...; kwargs...)

