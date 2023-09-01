# ------------------------------------------------------------------------
# API
# ------------------------------------------------------------------------
# Error for wrong number of input datasets.
test = SurrogateTest(MIShannon(), KSG1())
x, y, z = rand(30), rand(30), rand(30)
@test_throws ArgumentError independence(test, x)
@test_throws ArgumentError independence(test, x, y, z)


# Pairwise measures
include("MutualInformation.jl")
include("SMeasure.jl")
include("HMeasure.jl")
include("MMeasure.jl")
include("LMeasure.jl")
include("TransferEntropyPairwise.jl")
include("crossmappings.jl")

# Conditional measures
include("ConditionalMutualInformation.jl")
include("pmi.jl")

# Pairwise + conditional
include("TransferEntropyConditional.jl")
