"""
    SubsampledDataCausalityTest

A causality test where the data is subsampled in some manner before
applying the test.
"""
abstract type SubsampledDataCausalityTest{CT} <: MetaCausalityTest{CT} end 

get_ηs(x::SubsampledDataCausalityTest) = get_ηs(x.test)

export SubsampledDataCausalityTest, get_ηs

include("RandomSequencesTest.jl")
include("causality_RandomSequencesTest.jl")