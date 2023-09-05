# This file extends the single-variable information API in ComplexityMeasures.jl.

include("information_measures.jl")
include("information_functions.jl")

# Convenience encodings (should be in ComplexityMeasures.jl).
include("encoding/categorical_encoding.jl")

# Concrete ways of encoding multiple input datasets.
include("encoding/per_point_encoding.jl")
