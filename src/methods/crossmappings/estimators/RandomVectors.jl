using Random

export RandomVectors

"""
    RandomVectors <: CrossmapEstimator
    RandomVectors(; libsizes, replace = false, rng = Random.default_rng())

Cross-map over `N` different libraries, where `N = length(libsizes)`, and the `i`-th
library has cardinality `k = libsizes[i]`. Points within each library are randomly
selected, independently of other libraries, and `replace` controls whether or not to
sample with replacement. A user-specified `rng` may be specified for reproducibility.

This is method 3 from Luo et al. (2015)[^Luo2015].

[^Luo2015]:
    "Questionable causality: Cosmic rays to temperature." Proceedings of the National
    Academy of Sciences Aug 2015, 112 (34) E4638-E4639; DOI: 10.1073/pnas.1510571112
    Ming Luo, Holger Kantz, Ngar-Cheung Lau, Wenwen Huang, Yu Zhou.

See also: [`CrossmapEstimator`](@ref).
"""
struct RandomVectors{I, R} <: CrossmapEstimator{I, R}
    libsizes::I
    rng::R
    replace::Bool
    function RandomVectors(; libsizes::I, replace::Bool = false,
            rng::R = Random.default_rng()) where {I,R}
        new{I, R}(libsizes, rng, replace)
    end
end
