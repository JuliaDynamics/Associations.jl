include("core.jl")

# Counting for contingency matrices.
include("contingency_table.jl")

# Convenience encodings (should be in ComplexityMeasures.jl).
include("encoding/categorical_encoding.jl")

# Concrete ways of encoding multiple input datasets.
include("encoding/per_point_encoding.jl")
