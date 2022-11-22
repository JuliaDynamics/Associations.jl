using StaticArrays
using DelayEmbeddings

"""
    spikeembed(i, iₜmin, x; spike = one(eltype(x)), m::Int = 2) where T <: Number

Given scalar-valued `x`, and history length `m`, make a *continuous-time* embedding,
by counting the duration between spikes relative to time index `iₜ`.

`iₜmin` is the  minimum index for which to start checking, to ensure that a
length-`m` embedding vector is actually possible to construct. This must have been
checked beforehand by using [`first_index`](@ref), otherwise errors will be thrown.

See Fig 10 in Shorten et al. (2021) for details.

!!! note "Availability of length-`m` embedding vectors"
    This function assumes that there are at least `m` spikes available before time index
    `iₜ`. This function does not check for the existence of such spikes. Pre-check by using
    [`first_index`](@ref). If no spikes exist to form an embedding, an error will be thrown.

[^Shorten2021]: Shorten, D. P., Spinney, R. E., & Lizier, J. T. (2021). Estimating transfer
    entropy in continuous time between neural spike trains or other event-based data. PLoS
    computational biology, 17(4), e1008054.
"""
function spikeembed(target::AbstractVector{T}, event::EventIdentifier = OneSpike();
            m::Int = 2, # The same embedding for all timeseries processes for now.
            ) where T
    spike_identifier = one(T)

    # We can't know a priori how many elements this vector will contain,
    # because we're counting time *intervals* relative to some tᵢ, not
    # *values at particular times relative to tᵢ*.
    emb = Vector{SVector{m, T}}(undef, 0)
    iₜmin = first_index(target; m)

    # Pre-allocate single vector to avoid excessive allocations, and convert
    # back to SVector inside the inner function. If the compiler is feeling fine, then
    # this *could* incur no extra cost.
    p = zeros(MVector{m, T})

    L = length(target)
    for iₜ in iₜmin:L
        push!(emb, embed_at_iₜ_event!(p, iₜ, target, spike_identifier; m))
    end
    return emb
end

# Embeddings vectiors constructed at *all* times (independent of target spiking events).
function embed_at_alltimes(x::AbstractDataset{D, T}, event::EventIdentifier = OneSpike();
        m::Int = 2) where {D, T}
    spike_identifier = one(T)

    # Here, we don't care about the timing of spikes in the target time series,
    # so construct embedding points for all possible time indices `iₜ`
    iₜmin = first_index(x; m)
    L = length(x)

    # After finding the minimum required index for `m`-length history vectors to
    # exist, we can also know the size of the embeddings.
    embeddings = [Vector{SVector{m, Int}}(undef, (L - iₜmin + 1)) for d = 1:D]

    # Pre-allocate single vector to avoid excessive allocations, and convert
    # back to SVector inside the inner function. If the compiler is feeling fine, then
    # this *could* incur no extra cost.
    p = zeros(MVector{m, T})

    for d = 1:D
        ts = x[:, SVector{1, Int}(d)]
        for (i, iₜ) in enumerate(iₜmin:L)
            embeddings[d][i] = embed_at_iₜ_event!(p, iₜ, ts, spike_identifier; m)
        end
    end

    return Dataset.(embeddings)
end


# Embeddings constructed only *at* target spikes.
function embed_at_targetspikes(x::AbstractDataset{D, T}, event::EventIdentifier = OneSpike();
        m::Int = 2) where {D, T}
    spike_identifier = one(T)

    # Spikes in the *target* time series control which time indices for which
    # we construct embeddings.
    target = x[:, 1]
    inds_spikes = findall(target .== spike_identifier)

    # Ensure that all time series have enough spikes to be embedded.
    iₜmin = first_index(x; m)
    inds_spikes_valid = inds_spikes[inds_spikes .> iₜmin]
    L = length(inds_spikes_valid)

    # After finding the minimum required index for `m`-length history vectors to
    # exist, we can also know the size of the embeddings.
    embeddings = [Vector{SVector{m, Int}}(undef, L) for d = 1:D]

    # Pre-allocate single vector to avoid excessive allocations, and convert
    # back to SVector inside the inner function. If the compiler is feeling fine, then
    # this *could* incur no extra cost.
    p = zeros(MVector{m, T})

    for d in 1:D
        ts = x[:, SVector{1, Int}(d)]
        for (i, iₜ) in enumerate(inds_spikes_valid)
            embeddings[d][i] = embed_at_iₜ_event!(p, iₜ, ts, spike_identifier; m)
        end
    end
    return Dataset.(embeddings)
end

function spike_embeddings(x::AbstractDataset{D, T}, event::EventIdentifier = OneSpike();
        m::Int = 2) where {D, T}
    Eₜ = embed_at_targetspikes(x, event; m)
    Eᵤ = embed_at_alltimes(x, event; m)
    if D <= 2
        vars_cond = 1:2
    else
        vars_cond = [1; 3:length(Eₜ)]
    end
    # Joint embedding vector sets J include all variables, while conditional embedding sets
    # vector sets C excludes the source variable, which is assumed to always be
    # in position 2.
    Jₓ = hcat(Eₜ...)
    Cₓ = hcat(Eₜ[vars_cond]...)
    Jᵤ = hcat(Eᵤ...)
    Cᵤ = hcat(Eᵤ[vars_cond]...)

    return Jₓ, Cₓ, Jᵤ, Cᵤ
end

# Given a time index iₜ, find the `m` first intra-spike time intervals between
# t = t(iₜ) and t = t(1). Store the result in the pre-allocated MVector `v`.
# Spikes are identified using `spike_identifier` (usually `1`), but can
# be any identifier. Assumes `x` is a view of the [1:iₜ] first elements of x (Vector{StaticVector}).
function embed_at_iₜ_event!(v, iₜ, x, event_identifier; m = 2)
    n_found::Int = 0 # The number of embedding entries found.
    n_checked::Int = 0 # How many checked since last spike was found?
    k::Int = iₜ

    @inbounds while n_found < m
        if x[k - 1][] == event_identifier
            n_found += 1
            v[n_found] = n_checked
            n_checked = 0 # Restart counting from the current spike.
        else
            n_checked += 1
        end
        k -= 1
    end

    return SVector(v)
end



# using Test
# using DelayEmbeddings
# zt = [1, 0, 1, 0, 0, 1, 1, 0, 1, 0]
# xt = [0, 1, 0, 1, 0, 1, 0, 0, 0, 1]
# yt = [1, 1, 0, 0, 1, 0, 1, 0, 0, 1]
# Dt = Dataset(xt, yt, zt)

# # If we encounter two successive spikes, then that is counted as a 0 time interval.
# expected_m2 = [(0, 1), (1, 1), (2, 1), (0, 2), (0, 0), (1, 0), (0, 1)]
# expected_m3 = [(0, 2, 1), (0, 0, 2), (1, 0, 2), (0, 1, 0)]
# expected_m4 = [(0, 0, 2, 1), (1, 0, 2, 1), (0, 1, 0, 2)]
# @test spikeembed(zt, m = 2) == SVector{2, Int}.(expected_m2)
# @test spikeembed(zt, m = 3) == SVector{3, Int}.(expected_m3)
# @test spikeembed(zt, m = 4) == SVector{4, Int}.(expected_m4)

# # At all time_points
# @test embed_at_alltimes(Dt; m = 2) == Dataset.([
#     SVector{2, Int32}.([[0, 1], [1, 1], [0, 1], [1, 1], [2, 1], [3, 1]]),
#     SVector{2, Int32}.([[2, 0], [0, 2], [1, 2], [0, 1], [1, 1], [2, 1]]),
#     SVector{2, Int32}.([[1, 1], [2, 1], [0, 2], [0, 0], [1, 0], [0, 1]]),
# ])

# # At events for specific time series
# @test embed_at_targetspikes(Dt; m = 2) == Dataset.([
#     SVector{2, Int32}.([[1, 1], [3, 1]]),
#     SVector{2, Int32}.([[0, 2], [2, 1]]),
#     SVector{2, Int32}.([[2, 1], [0, 1]]),
# ])
