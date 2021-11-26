export pai

"""
    pai(x, y, d, τ; w = 0, correspondence_measure = Statistics.cor) → Float64
    pai(x, y, d, τ, bootstrap_method::Symbol; w = 0, correspondence_measure = Statistics.cor,
        method = :segment, L = ceil(Int, (length(x)-d*τ)*0.2), nreps = 100) → Vector{Float64}

Compute the pairwise asymmetric inference (PAI; McCracken, 2014)[^McCracken2014] between `x` and `y`.
Returns the correspondence between original and cross mapped values (the default is 
`ρ = correspondence_measure(y(t), ỹ(t) | M_xy)`).

PAI is a modification to Sugihara et al. (2012)'s CCM algorithm[^Sugihara2012], where instead of 
using completely out-of-sample prediction when trying to predict ``y(t)``, values about *both* variables 
are included in the embedding used to make predictions. Specifically, PAI computes the 
correspondence between the values ``y(t)`` and the cross-map estimated values ``ỹ(t) | M_xy``,
where the ``\\tilde{y}(t)`` are the values estimated using the embedding ``M_{xy} = \\{ ( x_t, x_{t-\\tau}, x_{t-2\\tau}, \\ldots, x_{t-(d - 1)\\tau} ) \\}``.
*Note: a `d+1`-dimensional embedding is used, rather than the `d`-dimensional embedding used for CCM.

Like for the CCM algorithm, the Theiler window `w` indicates how many temporal neighbors of the predictee is to be excluded 
during the nearest neighbors search (the default `w = 0` excludes only the predictee itself, while 
`w = 2` excludes the point itself plus its two nearest neighbors in time).

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
function pai(x, y, d, τ; correspondence_measure = Statistics.cor, w = 0) 
    Mₓy = crossmapembed(x, y, d, τ, PAIEmbedding()) 

    theiler = Theiler(w); 
    tree = KDTree(Mₓy)
    idxs = bulkisearch(tree, Mₓy, NeighborNumber(d+1), theiler)
    ỹ = copy(y)

    for i in 1:length(Mₓy)
        J = idxs[i]
        xᵢ = Mₓy[i]
        n1 = norm(xᵢ - Mₓy[J[1]])
        # If distance to nearest neighor is zero, then we get division by zero.
        # If that is the case, set all weights all to zero.
        if n1 == 0.0
            w .= 0.0
        else
            w = [exp(-norm(xᵢ - Mₓy[j])/n1) for j in J]
            w ./= sum(w)
        end
        ỹ[i] = sum(w[k]*y[j] for (k, j) in enumerate(J))
    end

    return correspondence_measure(y, ỹ)
end


function pai(x, y, d, τ, bootstrap_method::Symbol;
        L = ceil(Int, (length(x) - d * τ) * 0.2), nreps = 100, 
        w = 0, correspondence_measure = Statistics.cor)
    
    @assert length(x) == length(y)
    @assert length(x) - d * τ > 0

    Mₓy = crossmapembed(x, y, d, τ, PAIEmbedding()) 
    theiler = Theiler(w);

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
            idxs = StatsBase.sample(1:length(Mₓy), L)
        elseif bootstrap_method == :segment 
            # Training set is a random set of L embedding vectors in Mₓ
            # This is method 2 in Luo et al. (2015), but with added exclusion of 
            # the predictee from the libraries.
            startidx = rand(1:length(Mₓy) - L)
            idxs = startidx:startidx + L - 1
        else
            error("$bootstrap_method is not a valid bootstrap method")
        end
        training_set = Mₓy[idxs]
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
            # division by zero when computing weights. In that case, we set all weights to zero.
            dist₁ = norm(xᵢ - nns[1])

            if dist₁ == 0.0
                w .= 0.0
            else 
                for j = 1:d+1
                    distⱼ = norm(xᵢ - nns[j])
                    u[j] = exp(-distⱼ / dist₁)
                end
                for j = 1:d+1
                    w[j] = u[j] / sum(u)
                end
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