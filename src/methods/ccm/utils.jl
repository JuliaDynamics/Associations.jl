using Neighborhood: bulkisearch
using Neighborhood: NeighborNumber, Theiler, KDTree
using DelayEmbeddings: Dataset
using DelayEmbeddings: genembed
using Statistics: cor
using StatsBase: sample
using LinearAlgebra: norm

abstract type CrossmapEmbedding end

struct CCMEmbedding end
struct PAIEmbedding end

function crossmapembed(x, d, τ, method::CCMEmbedding)
    τs = 0:-τ:d*-τ+1
    Mₓ = genembed(x, τs);
    return Mₓ
end

function crossmapembed(x, y, d, τ, method::PAIEmbedding)
    τs = [collect(0:-τ:d*-τ+1); 0]
    js = [repeat([1], d); 2]
    Mₓy = genembed(Dataset(x, y), τs, js);
    return Mₓy
end

function crossmap_basic(Mₓ, y, d, τ; correspondence_measure = cor, r = 0)
    theiler = Theiler(r);
    tree = KDTree(Mₓ)
    idxs = bulkisearch(tree, Mₓ, NeighborNumber(d+1), theiler)
    ỹ = copy(y)

    for i in eachindex(Mₓ)
        J = idxs[i]
        xᵢ = Mₓ[i]
        n1 = norm(xᵢ - Mₓ[J[1]])

        # If distance to nearest neighor is zero, then we get division by zero.
        # If that is the case, set all weights all to zero (or just skip
        # computing ỹᵢ, because it is already set to 0).
        n1 > 0.0 || continue
        w = [exp(-norm(xᵢ - Mₓ[j])/n1) for j in J]
        w ./= sum(w)
        ỹ[i] = sum(w[k]*y[j] for (k, j) in enumerate(J))
    end

    return correspondence_measure(y, ỹ)
end

function crossmap_bootstrap(Mₓ, y, d, τ, bootstrap_method::Symbol;
    L = ceil(Int, (length(y) - d * τ) * 0.2), nreps = 100,
    r = 0, correspondence_measure = cor)

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
            idxs = sample(1:length(Mₓ), L)
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

        @inbounds for i in eachindex(training_set)
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
            s = sum(u)
            for j = 1:d+1
                w[j] = u[j] / s
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
