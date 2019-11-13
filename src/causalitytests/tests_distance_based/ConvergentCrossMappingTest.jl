
"""
    ConvergentCrossMappingTest(timeseries_lengths; 
        dim::Int = 3, τ::Int = 1, libsize::Int = 10,
        replace::Bool = false, n_reps::Int = 100, surr_func::Function = randomshuffle,
        which_is_surr::Symbol = :none, exclusion_radius::Int = 0,
        tree_type = NearestNeighbors.KDTree, distance_metric = Distances.Euclidean(),
        correspondence_measure = StatsBase.cor, 
        η::Int = 0)

The parameters for a convergent cross mapping [1] test.

## Mandatory keyword arguments 

- **`timeseries_lengths`**: The time series lengths over which to cross map and check convergence.

## Optional keyword arguments 

- **`dim`**: The dimension of the state space reconstruction (delay embedding)
    constructed from the `response` series. Default is `dim = 3`.
- **`τ`**: The embedding lag for the delay embedding constructed from `response`.
    Default is `τ = 1`.
- **`η`**: The prediction lag to use when predicting scalar values of `driver`
    fromthe delay embedding of `response`.
    `η > 0` are forward lags (causal; `driver`'s past influences `response`'s future),
    and `η < 0` are backwards lags (non-causal; `driver`'s' future influences
    `response`'s past). Adjust the prediction lag if you
    want to performed lagged ccm
    [(Ye et al., 2015)](https://www.nature.com/articles/srep14750).
    Default is `η = 0`, as in
    [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079).
    *Note: The sign of the lag `η` is organized to conform with the conventions in
    [TransferEntropy.jl](), and is opposite to the convention used in the
    [`rEDM`](https://cran.r-project.org/web/packages/rEDM/index.html) package
    ([Ye et al., 2016](https://cran.r-project.org/web/packages/rEDM/index.html)).*
- **`libsize`**: Among how many delay embedding points should we sample time indices
    and look for nearest neighbours at each cross mapping realization (of which there
    are `n_reps`)?
- **`n_reps`**: The number of times we draw a library of `libsize` points from the
    delay embedding of `response` and try to predict `driver` values. Equivalently,
    how many times do we cross map for this value of `libsize`?
    Default is `n_reps = 100`.
- **`replace`**: Sample delay embedding points with replacement? Default is `replace = true`.
- **`exclusion_radius`**: How many temporal neighbors of the delay embedding
    point `response_embedding(t)` to exclude when searching for neighbors to
    determine weights for predicting the scalar point `driver(t + η)`.
    Default is `exclusion_radius = 0`.
- **`which_is_surr`**: Which data series should be replaced by a surrogate
    realization of the type given by `surr_type`? Must be one of the
    following: `:response`, `:driver`, `:none`, `:both`.
    Default is `:none`.
- **`surr_func`**: A valid surrogate function from TimeseriesSurrogates.jl.
- **`tree_type`**: The type of tree to build when looking for nearest neighbors.
    Must be a tree type from NearestNeighbors.jl. For now, this is either
    `BruteTree`, `KDTree` or `BallTree`.
- **`distance_metric`**: An instance of a `Metric` from Distances.jl. `BallTree` and `BruteTree` work with any `Metric`.
    `KDTree` only works with the axis aligned metrics `Euclidean`, `Chebyshev`,
    `Minkowski` and `Cityblock`. Default is `metric = Euclidean()` *(note the instantiation of the metric)*.
- **`correspondence_measure`**: The function that computes the correspondence
    between actual values of `driver` and predicted values. Can be any
    function returning a similarity measure between two vectors of values.
    Default is `correspondence_measure = StatsBase.cor`, which returns values on ``[-1, 1]``.
    In this case, any negative values are usually filtered out (interpreted as zero coupling) and
    a value of ``1`` means perfect prediction.
    [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079)
    also proposes to use the root mean square deviation, for which a value of ``0`` would
    be perfect prediction.
- **`summarise`**: Should cross map skills be summarised for each time series length? Default is 
    `summarise = false`. 
- **`average_measure`**: Either `:median` or `:mean`. Default is `:median`
- **`uncertainty_measure`**: Either `:quantile` or `:std`.
- **`quantiles`**: Compute uncertainty over quantile(s) if `uncertainty_measure` is `:quantile`. 
    Default is `[0.327, 0.673]`, roughly corresponding to `1s` for normally distributed data.

## References 

1. Sugihara, George, et al. "Detecting causality in complex ecosystems." Science (2012): 1227079. 
    [http://science.sciencemag.org/content/early/2012/09/19/science.1227079](http://science.sciencemag.org/content/early/2012/09/19/science.1227079)

2. Ye, Hao, et al. "Distinguishing time-delayed causal interactions using convergent cross mapping." Scientific 
    Reports 5 (2015): 14750. [https://www.nature.com/articles/srep14750](https://www.nature.com/articles/srep14750)

3. Ye, H., et al. "rEDM: Applications of empirical dynamic modeling from time series." R Package Version 
    0.4 7 (2016). [https://cran.r-project.org/web/packages/rEDM/index.html](https://cran.r-project.org/web/packages/rEDM/index.html)
"""
Base.@kwdef struct ConvergentCrossMappingTest{N, NL} <: DistanceBasedCausalityTest{N}
    """ 
    The dimension of the state space reconstruction (delay embedding) constructed 
    from the response series. Default is `dim = 3`. 
    """
    dim::Int = 3
    
    """ The embedding lag for the delay embedding constructed from response. Default is τ = 1. """
    τ::Int = 1
    
    """ 
    The prediction lag to use when predicting scalar values of `driver`
    from the delay embedding of `response`.
    `η > 0` are forward lags (causal; `driver`'s past influences `response`'s future),
    and `η < 0` are backwards lags (non-causal; `driver`'s' future influences
    `response`'s past). Adjust the prediction lag if you
    want to performed lagged ccm
    [(Ye et al., 2015)](https://www.nature.com/articles/srep14750).
    Default is `η = 0`, as in
    [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079).
    *Note: The sign of the lag `η` is organized to conform with the conventions in
    [TransferEntropy.jl](), and is opposite to the convention used in the
    [`rEDM`](https://cran.r-project.org/web/packages/rEDM/index.html) package
    ([Ye et al., 2016](https://cran.r-project.org/web/packages/rEDM/index.html)).*
    """ 
    η::Int = 0

    """ 
    Among how many delay embedding points should we sample time indices
    and look for nearest neighbours at each cross mapping realization (of which there
    are `n_reps`)? 
    """
    libsize::Int = 10
    
    """ Sample delay embedding points with replacement? Default is `replace = true`. """
    replace::Bool = true
    
    """
    The number of times we draw a library of `libsize` points from the
    delay embedding of `response` and try to predict `driver` values. Equivalently,
    how many times do we cross map for this value of `libsize`?
    Default is `n_reps = 100`.
    """
    n_reps::Int = 100
    
    """ A valid surrogate function from TimeseriesSurrogates.jl. """
    surr_func::Function = randomshuffle
    
    """ 
    Which data series should be replaced by a surrogate
    realization of the type given by `surr_type`? Must be one of the
    following: `:response`, `:driver`, `:none`, `:both`.
    Default is `:none`.
    """
    which_is_surr::Symbol = :none
    
    """ 
    How many temporal neighbors of the delay embedding
    point `response_embedding(t)` to exclude when searching for neighbors to
    determine weights for predicting the scalar point `driver(t + η)`.
    Default is `exclusion_radius = 0`.
    """
    exclusion_radius::Int = 0
    
    """
    The type of tree to build when looking for nearest neighbors.
    Must be a tree type from NearestNeighbors.jl. For now, this is either
    `BruteTree`, `KDTree` or `BallTree`.
    """
    tree_type = NearestNeighbors.KDTree
    
    """
    An instance of a `Metric` from Distances.jl. `BallTree` and `BruteTree` work with any `Metric`.
    `KDTree` only works with the axis aligned metrics `Euclidean`, `Chebyshev`,
    `Minkowski` and `Cityblock`. Default is `metric = Euclidean()` 
    *(note the instantiation of the metric)*.
    """
    distance_metric = Distances.Euclidean()

    """
    The function that computes the correspondence between actual values of `driver` and predicted 
    values. Can be any function returning a similarity measure between two vectors of values.
    Default is `correspondence_measure = StatsBase.cor`, which returns values on ``[-1, 1]``.
    In this case, any negative values are usually filtered out (interpreted as zero coupling) and
    a value of ``1`` means perfect prediction.
    [Sugihara et al. (2012)](http://science.sciencemag.org/content/early/2012/09/19/science.1227079)
    also proposes to use the root mean square deviation, for which a value of ``0`` would
    be perfect prediction.
    """
    correspondence_measure = StatsBase.cor

    """ Should cross map skills be summarised for each time series length? Default is `summarise = false`. """
    summarise::Bool = false 

    """ Either `:median` or `:mean`. Default is `:median`. """
    average_measure = :median

    """ Either `:quantile` or `:std`. """
    uncertainty_measure = :quantile

    """ 
    Compute uncertainty over quantile(s) if `uncertainty_measure` is `:quantile`. 
    Default is ``, roughly corresponding to `1s` for normally distributed data. 
    """
    quantiles = [0.327, 0.673] 

    """ The time series lengths over which to cross map and check convergence """
    timeseries_lengths

    function ConvergentCrossMappingTest(dim::Int, τ::Int, η::Int, libsize::Int, replace::Bool, 
        n_reps::Int, surr_func::Function, which_is_surr::Symbol, 
        exclusion_radius::Int, tree_type, distance_metric, correspondence_measure,
        summarise::Bool, average_measure, uncertainty_measure, quantiles, 
        timeseries_lengths)
        
        # The number of returned elements. We're performing the test for `n_reps` 
        # different randomly selected point libraries, so `N = n_reps`.
        N = n_reps

        # The number of different time series lengths 
        NL = length(timeseries_lengths)

        new{N, NL}(dim, τ, η, libsize, replace, n_reps, surr_func, which_is_surr, 
            exclusion_radius, tree_type, distance_metric, correspondence_measure,
            summarise, average_measure, uncertainty_measure, quantiles, 
            timeseries_lengths)
    end
end


function causality(source::AbstractVector{T}, target::AbstractVector{T}, p::ConvergentCrossMappingTest) where {T<:Real}
     convergentcrossmap(source, target, p.timeseries_lengths;
        dim = p.dim,
        τ = p.τ,
        ν = p.η,
        libsize = p.libsize,
        n_reps = p.n_reps,
        replace = p.replace,
        exclusion_radius = p.exclusion_radius,
        jitter = 1e-5 * max(StatsBase.std(source), StatsBase.std(target)),
        which_is_surr = p.which_is_surr,
        surr_func = p.surr_func,
        tree_type = p.tree_type,
        distance_metric = p.distance_metric,
        correspondence_measure = p.correspondence_measure,
        summarise = p.summarise,
        quantiles = p.quantiles,
        uncertainty_measure = p.uncertainty_measure,
        average_measure = p.average_measure)
end


export ConvergentCrossMappingTest