# Generic interface
include("te_from_marginal_entropy.jl")

# Specialized methods
include("ordinal/symbolic_te.jl")
include("hilbert/hilbert.jl")

# Advanced methods (automatic variable selection)
include("advanced/bbnue.jl")
