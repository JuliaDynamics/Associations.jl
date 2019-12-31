"""
    PreprocessedDataCausalityTest <: MetaCausalityTest

A causality test where the data is subjected to some preprocessing before 
applying the test.
"""
abstract type PreprocessedDataCausalityTest{CT} <: MetaCausalityTest{CT} end 

export PreprocessedDataCausalityTest

include("BinnedDataCausalityTest.jl")
include("ConstrainedTest.jl")
include("InterpolateBinTest.jl")

include("causality_BinnedDataCausalityTest.jl")
include("causality_ConstrainedTest.jl")
include("causality_InterpolateBinTest.jl")