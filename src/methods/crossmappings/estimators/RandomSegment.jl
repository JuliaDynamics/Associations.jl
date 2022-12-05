
"""
    RandomSegment <: CrossmapEstimator
    RandomSegment(L::Int)

Indicatates that cross mapping is performed on contiguous time series
segments/windows of length `L` with a randomly selected starting point.
This is method 2 from Luo et al. (2015)[^Luo2015].

[^Luo2015]:
    "Questionable causality: Cosmic rays to temperature." Proceedings of the National
    Academy of Sciences Aug 2015, 112 (34) E4638-E4639; DOI: 10.1073/pnas.1510571112
    Ming Luo, Holger Kantz, Ngar-Cheung Lau, Wenwen Huang, Yu Zhou.
"""
struct RandomSegment <: CrossmapEstimator
    L::Int
end
