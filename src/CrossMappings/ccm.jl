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
    Mₓ = embed(x, d, τ); 
    theiler = Theiler(r); 
    tree = KDTree(Mₓ)
    idxs = bulkisearch(tree, Mₓ, NeighborNumber(d+1), theiler)
    ỹ = zeros(length(y))
    
    for i in 1:length(Mₓ)
        J = idxs[i]
        xᵢ = Mₓ[i]
        n1 = norm(xᵢ - Mₓ[J[1]])
   
        # If distance to nearest neighor is zero, then we get division by zero.
        # If that is the case, set all weights all to zero (or just skip 
        # computing ỹᵢ, because it is already set to 0).
        n1 > 0.0 || continue
        w = [exp(-norm(xᵢ - Mₓ[j])/n1) for j in J]
        w ./= sum(w)
        ỹ[i] = sum(w[k]*y[j] for (k, j) in enumerate(J))
    end
    
    return correspondence_measure(y, ỹ)
end

function crossmap(x, y, d, τ, bootstrap_method::Symbol;
        L = ceil(Int, (length(x) - d * τ) * 0.2), nreps = 100, 
        r = 0, correspondence_measure = Statistics.cor)
    
    @assert length(x) == length(y)
    @assert length(x) - d * τ > 0

    Mₓ = embed(x, d, τ); 
    theiler = Theiler(r);

    # Compute correlations between out-of-sample target and embedding for `nreps` 
    # different training sets
    cors = zeros(nreps)

    # These arrays are re-used.
    u = zeros(d + 1)
    w = zeros(d + 1)
    ỹs = zeros(L)

    for n = 1:nreps
        # Select training set and keep track of corresponding y-values
        if bootstrap_method == :random 
            # Training set is a random embedding vectors. This is method 3 in 
            # Luo et al. (2015)
            idxs = StatsBase.sample(1:length(Mₓ), L)
        elseif bootstrap_method == :segment 
            # Training set is a random set of L embedding vectors in Mₓ
            # This is method 2 in Luo et al. (2015), but with added exclusion of 
            # the predictee from the libraries.
            startidx = rand(1:length(Mₓ) - L)
            idxs = startidx:startidx + L - 1
        else
            error("$bootstrap_method is not a valid bootstrap method")
        end
        training_set = Mₓ[idxs]
        ys = y[idxs]

        # Find the nearest neighbors of points in the training set
        tree = KDTree(training_set)
        idxs_nns = bulkisearch(tree, training_set, NeighborNumber(d+1), theiler)
        
        # Reset predictions
        ỹs .= 0.0

        @inbounds for i in 1:length(training_set) 
            # For every point xᵢ ∈ training_set
            xᵢ = training_set[i]

            # Get the indices of its nearest neighbors (in the training set, not the entire Mx)
            nns_idxs = idxs_nns[i]
            nns = training_set[nns_idxs]

            # Depending on the input data, we might encounter a case where the distance from xᵢ
            # to its closest neighbor is indistinguishable from zero. In that case, we encounter
            # division by zero when computing weights. In that case, set weights to zero (but in
            # practice, just continue to the next step, because ỹᵢ will is already initialized to 0)
            dist₁ = norm(xᵢ - nns[1])
            dist₁ > 0.0 || continue

            for j = 1:d+1
                distⱼ = norm(xᵢ - nns[j])
                u[j] = exp(-distⱼ / dist₁)
            end
            for j = 1:d+1
                w[j] = u[j] / sum(u)
            end
            
            # For each weight
            for j = 1:d+1
                # Find the scalar value corresponding to this weight
                idx_nnⱼ = nns_idxs[j]
                y_scalar = ys[idx_nnⱼ]
                ỹs[i] += w[j] * y_scalar
            end
        end
        cors[n] = correspondence_measure(ys, ỹs)
    end

    return cors
end