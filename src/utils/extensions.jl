import StateSpaceSets: AbstractDataset

Base.similar(x::T) where T <: AbstractDataset = T(similar(x.data))


# TODO: implement fast sampling for Datasets. There's a bit of work involved


# import StatsBase: direct_sample!, sample!
# using StatsBase: Sampler
# export direct_sample!
# function direct_sample!(rng::AbstractRNG, a::AbstractDataset, x::AbstractDataset)
#     s = Sampler(rng, 1:length(a))
#     for i in eachindex(x)
#         @inbounds x[i] = a[rand(rng, s)]
#     end
#     return x
# end
# direct_sample!(a::AbstractDataset, x::AbstractDataset) = direct_sample!(Random.GLOBAL_RNG, a, x)

# We're not returning ordered samples, so no need to override `ordered_sample!`

# # using Random: AbstractRNG
# function sample!(rng::AbstractRNG, a::AbstractDataset, x::AbstractDataset, args...; kwargs...)
#     sample!(rng, )
# end
