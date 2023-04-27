import StateSpaceSets: AbstractStateSpaceSet

Base.similar(x::T) where T <: AbstractStateSpaceSet = T(similar(x.data))


# TODO: implement fast sampling for StateSpaceSets. There's a bit of work involved


# import StatsBase: direct_sample!, sample!
# using StatsBase: Sampler
# export direct_sample!
# function direct_sample!(rng::AbstractRNG, a::AbstractStateSpaceSet, x::AbstractStateSpaceSet)
#     s = Sampler(rng, 1:length(a))
#     for i in eachindex(x)
#         @inbounds x[i] = a[rand(rng, s)]
#     end
#     return x
# end
# direct_sample!(a::AbstractStateSpaceSet, x::AbstractStateSpaceSet) = direct_sample!(Random.GLOBAL_RNG, a, x)

# We're not returning ordered samples, so no need to override `ordered_sample!`

# # using Random: AbstractRNG
# function sample!(rng::AbstractRNG, a::AbstractStateSpaceSet, x::AbstractStateSpaceSet, args...; kwargs...)
#     sample!(rng, )
# end
