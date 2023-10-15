using Random

export RandomVectors

"""
    RandomVectors <: CrossmapEstimator
    RandomVectors(definition::CrossmapMeasure; 
        libsizes, replace = false, rng = Random.default_rng())

Cross-map over `N` different libraries, where `N = length(libsizes)`, and the `i`-th
library has cardinality `k = libsizes[i]`. Points within each library are randomly
selected, independently of other libraries, and `replace` controls whether or not to
sample with replacement. A user-specified `rng` may be specified for reproducibility.

This is method 3 from [Luo2015](@citet).

See also: [`CrossmapEstimator`](@ref).
"""
struct RandomVectors{M <: CrossmapMeasure, I, R} <: CrossmapEstimator{M, I, R}
    definition::M
    libsizes::I
    rng::R
    replace::Bool
    function RandomVectors(definition::M; libsizes::I, replace::Bool = false,
            rng::R = Random.default_rng()) where {M<:CrossmapMeasure, I, R}
        new{M, I, R}(definition, libsizes, rng, replace)
    end
end

function library_indices(est::RandomVectors, i::Int, target, args...)
    N = length(target)
    L = est.libsizes[i]
    return library_indices(measure, est, N, L)
end

function library_indices(est::RandomVectors, N::Int, L::Int)
    return sample(est.rng, 1:N, L; replace = est.replace)
end
