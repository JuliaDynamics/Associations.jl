using CausalityTools,StatsBase, DelayEmbeddings, LinearAlgebra, Statistics
import DelayEmbeddings.Neighborhood: Theiler, bulkisearch, NeighborNumber

export crossmap

"""
    crossmap(x, y, d, τ; r = 0, correspondence_measure = Statistics.cor) → Float64
    crossmap(x, y, d, τ, bootstrap_method::Symbol; r = 0, correspondence_measure = Statistics.cor,
        method = :segment, L = ceil(Int, (length(x)-d*τ)*0.2), nreps = 100) → Vector{Float64}

Compute the cross mapping [^Sugihara2012] between `x` and `y`, which is the correspondence (computed using 
`correspondence measure`) between the values ``y(t)`` and the cross-map estimated values ``ỹ(t) | M_x``.
Returns the correspondence between original and cross mapped values (the default is 
`ρ = correspondence_measure(y(t), ỹ(t) | M_x)`).

Here, ``y(t)`` are the raw values of the time series `y`, and ``ỹ(t)`` are the predicted values 
computed from the out-of-sample embedding ``M_X`` constructed from the time series `x` with 
embedding dimension `d` and embedding lag `τ`.

The Theiler window `r` indicates how many temporal neighbors of the predictee is to be excluded 
during the nearest neighbors search (the default `r = 0` excludes only the predictee itself, while 
`r = 2` excludes the point itself plus its two nearest neighbors in time).

If `bootstrap_method` is specified, then `nreps` different bootstrapped estimates of 
`correspondence_measure(y(t), ỹ(t) | M_x)` are returned. The following bootstrap methods are available:

- `bootstrap_method = :random` selects training sets of length `L` consisting of randomly selected 
    points from the embedding ``M_x``  (time ordering does not matter). This is method 3 from Luo 
    et al. (2015)[^Luo2015], which critiqued the original Sugihara et al. methodology.
- `bootstrap_method = :segment` selects training sets consisting of time-contiguous segments 
    (each of lenght `L`) of embedding vectors in ``M_x`` (time ordering matters). This is 
    method 2 from Luo et al. (2015)[^Luo2015].

[^Sugihara2012]: Sugihara, George, et al. "Detecting causality in complex ecosystems." Science (2012): 1227079.[http://science.sciencemag.org/content/early/2012/09/19/science.1227079](http://science.sciencemag.org/content/early/2012/09/19/science.1227079)
[^Luo2015]: "Questionable causality: Cosmic rays to temperature." Proceedings of the National Academy of Sciences Aug 2015, 112 (34) E4638-E4639; DOI: 10.1073/pnas.1510571112 Ming Luo, Holger Kantz, Ngar-Cheung Lau, Wenwen Huang, Yu Zhou
"""
function crossmap(x, y, d, τ; correspondence_measure = Statistics.cor, r = 0) 
    # Embed `x`
    Mₓ = crossmapembed(x, d, τ, CCMEmbedding())

    # Compute `ỹ | Mₓ` and return `correspondence_measure(ỹ | Mₓ, y).
    return crossmap_basic(Mₓ, y, d, τ, correspondence_measure = correspondence_measure, r = r)
end

function crossmap(x, y, d, τ, bootstrap_method::Symbol;
        L = ceil(Int, (length(x) - d * τ) * 0.2), nreps = 100, 
        r = 0, correspondence_measure = Statistics.cor)
    
    # Embed `x`
    @assert length(x) == length(y)
    @assert length(x) - d * τ > 0
    Mₓ = crossmapembed(x, d, τ, CCMEmbedding());

    # Compute `ỹ | Mₓ` and return `nreps` values for `correspondence_measure(ỹ | Mₓ, y),
    # each value being computed on a random library (sequential or not) consisting of `L` points.
    return crossmap_bootstrap(Mₓ, y, d, τ, bootstrap_method; L = L, nreps = nreps, r = r, 
        correspondence_measure = correspondence_measure)
end