using Statistics
export pai

"""
    crossmap(x, y, d, τ; r = 0, correspondence_measure = Statistics.cor) → Float64
    crossmap(x, y, d, τ, bootstrap_method::Symbol; r = 0, correspondence_measure = Statistics.cor,
        method = :segment, L = ceil(Int, (length(x)-d*τ)*0.2), nreps = 100) → Vector{Float64}

!!! info "This syntax is deprecated"
    This syntax is deprecated. It will continue to work for CausalityTools v1.X, but will
    be removed in CausalityTools v2. See [here](@ref cross_mapping_api) for updated syntax.

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
    @warn """\
    `crossmap(x, y, d, τ; kwargs...)`d is deprecated. Use \
    `crossmap(CCM(; d, τ), x, y` instead.
    """

    measure = CCM(; d, τ, f = correspondence_measure, w = r)
    return crossmap(measure, y, x)
end

function crossmap(x, y, d, τ, bootstrap_method::Symbol;
        L = ceil(Int, (length(x) - d * τ) * 0.2), nreps = 100,
        r = 0, correspondence_measure = Statistics.cor)
    @warn """\
    `crossmap(x, y, d, τ, bootstrap_method; kwargs...)` is deprecated. Use \
    `crossmap(Ensemble(PAI(; d, τ), est::CrossmapEstimator), x, y)` instead.
    """
    measure = ConvergentCrossMapping(; d, τ, f = correspondence_measure, w = r)
    if bootstrap_method == :random
        est = RandomVectors(libsizes = L)
    else bootstrap_method == :segment
        est = RandomSegment(libsizes = L)
    end
    ensemble = Ensemble(measure, est, nreps = nreps)
    return crossmap(ensemble, y, x)
end

"""
    pai(x, y, d, τ; w = 0, correspondence_measure = Statistics.cor) → Float64
    pai(x, y, d, τ, bootstrap_method::Symbol; w = 0, correspondence_measure = Statistics.cor,
        method = :segment, L = ceil(Int, (length(x)-d*τ)*0.2), nreps = 100) → Vector{Float64}

!!! info "This syntax is deprecated"
    This syntax is deprecated. It will continue to work for CausalityTools v1.X, but will
    be removed in CausalityTools v2. See [here](@ref cross_mapping_api) for updated syntax.

Compute the pairwise asymmetric inference (PAI; McCracken, 2014)[^McCracken2014] between `x` and `y`.
Returns the correspondence between original and cross mapped values (the default is
`ρ = correspondence_measure(y(t), ỹ(t) | M_xy)`).

PAI is a modification to Sugihara et al. (2012)'s CCM algorithm[^Sugihara2012], where instead of
using completely out-of-sample prediction when trying to predict ``y(t)``, values about *both* variables
are included in the embedding used to make predictions. Specifically, PAI computes the
correspondence between the values ``y(t)`` and the cross-map estimated values ``ỹ(t) | M_xy``,
where the ``\\tilde{y}(t)`` are the values estimated using the embedding ``M_{xy} = \\{ ( x_t, x_{t-\\tau}, x_{t-2\\tau}, \\ldots, x_{t-(d - 1)\\tau} ) \\}``.
*Note: a `d+1`-dimensional embedding is used, rather than the `d`-dimensional embedding used for CCM.
Like for the CCM algorithm, the Theiler window `r` indicates how many temporal neighbors of the predictee is to be excluded
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

[^McCracken2014]: McCracken, James M., and Robert S. Weigel. "Convergent cross-mapping and pairwise asymmetric inference." Physical Review E 90.6 (2014): 062903.
[^Sugihara2012]: Sugihara, George, et al. "Detecting causality in complex ecosystems." Science (2012): 1227079.[http://science.sciencemag.org/content/early/2012/09/19/science.1227079](http://science.sciencemag.org/content/early/2012/09/19/science.1227079)
[^Luo2015]: "Questionable causality: Cosmic rays to temperature." Proceedings of the National Academy of Sciences Aug 2015, 112 (34) E4638-E4639; DOI: 10.1073/pnas.1510571112 Ming Luo, Holger Kantz, Ngar-Cheung Lau, Wenwen Huang, Yu Zhou
"""
function pai(x, y, d, τ; correspondence_measure = Statistics.cor, r = 0)
    @warn """\
    `pai(x, y, d, τ; kwargs...)`d is deprecated. Use \
    `crossmap(PAI(; d, τ), x, y)` instead.
    """

    measure = PairwiseAsymmetricInference(; d, τ, f = correspondence_measure, w = r)
    return crossmap(measure, y, x)
end

function pai(x, y, d, τ, bootstrap_method::Symbol;
        L = ceil(Int, (length(x) - d * τ) * 0.2), nreps = 100,
        r = 0, correspondence_measure = Statistics.cor)
    @warn """\
    `pai(x, y, d, τ; kwargs...)`d is deprecated. Use \
    `crossmap(Ensemble(PAI(; d, τ), est::CrossmapEstimator), x, y)` instead.
    """
    measure = CCM(; d, τ, f = correspondence_measure, w = r)
    if bootstrap_method == :random
        est = RandomVectors(libsizes = L)
    else bootstrap_method == :segment
        est = RandomSegment(libsizes = L)
    end
    ensemble = Ensemble(measure, est, nreps = nreps)
    return crossmap(ensemble, y, x)
end
