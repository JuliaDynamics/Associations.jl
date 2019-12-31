export MetaCausalityTest

abstract type MetaCausalityTest{CT <: CausalityTest} end

include("PreprocessedDataCausalityTest/PreprocessedDataCausalityTest.jl")
include("SubsampledDataCausalityTest/SubsampledDataCausalityTest.jl")



