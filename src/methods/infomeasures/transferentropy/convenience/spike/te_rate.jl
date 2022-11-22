
function event_transferentropyrate(e::Entropy,
        event::EventIdentifier,
        est::ProbabilitiesEstimator,
        target::V, source::V, cond::Vararg{V};
        m = 2) where V <: AbstractVector

    D = Dataset(target, source, cond...)

    # Construct embeddings, given that the events are identified by `event`.
    Jₓ, Cₓ, Jᵤ, Cᵤ = spike_embeddings(D, event; m)
    return entropy(e, Cₓ, est) +
        entropy(e, Jₓ, est) +
        entropy(e, Jᵤ, est) -
        entropy(e, Cᵤ, est)
end

transferentropy_event(event::EventIdentifier, est::ProbabilitiesEstimator,
    target::V, source::V, cond::Vararg{V};
    m = 2, base = 2) where {V <: AbstractVector} =
        transferentropy_event(Shannon(; base), event, est, target, source, cond; m)

est = CountOccurrences()
