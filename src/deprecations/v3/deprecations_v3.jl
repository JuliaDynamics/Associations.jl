# Cross mappings
include("crossmappings.jl")

# Closeness
include("joint_distance_distribution.jl")
include("smeasure.jl")
include("hmeasure.jl")
include("mmeasure.jl")
include("lmeasure.jl")

# Correlation
include("pearson_correlation.jl")
include("distance_correlation.jl")
include("partial_correlation.jl")

# Full names
@deprecate PMI PartialMutualInformation
