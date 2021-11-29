using Distances, Statistics, LinearAlgebra, DelayEmbeddings, Neighborhood, Distances, LinearAlgebra

export smap

"""
    smap(X::Dataset; θ = 1.0, k = 1, 
        training = 1:length(x) ÷ 2 - k, 
        predictees = length(x) ÷ 2 + 1:length(x)) → Vector{Float64}, Vector{Float64}

Sequential locally weighted global linear map (S-map; Sugihara, 1994)[^Sugihara1994].

For each predictee ``x_i ∈ X_{pred}```, the algorithm uses points in the library/training set 
``X_{train} \\setminus x_i`` to fit a weighted linear model for predicting ``x_{i+k}`` 
(xᵢ projected `k` time steps into the future).

Returns two scalar vectors: `x̂s = [x̂₁₊ₖ, x̂₂₊ₖ, ..., x̂ₙ+k]`, which are the predicted values, and 
`x̂s_truths = [x₁₊ₖ, x₂₊ₖ, ..., xₙ+k]`, which are the actual values ``x_{i+k} \\in X_{pred}`` 
for ``i = 1, \\ldots, n``, where ``n`` is the number of predictees.

## Input data

The algorithm approximates a global map based on *embedding points* and predicts scalar values 
in the *first* column of `X`. If your input data is a scalar time series, it must therefore 
first be embedded using `DynamicalSystems.genembed` first.

## Training and prediction sets

By default, the first `ntrain` points of `x` is used as the library set (`Xtrain = x[1:Ntrain-k]`), and the 
remaining half (`Xpred = x[Ntrain+1:end]`) is assigned to the prediction set. Overlapping index ranges 
are not possible as of yet.

When `θ = 0.0`, all weights are identical, and the A reduces to a linear autoregressive A. 
Nonlinearity is introduced when `θ > 0`, so tuning this parameters can be used to distinguish 
nonlinear dynamical systems from linear stochastic systems.

[^Sugihara1994]: Sugihara, G. (1994). Nonlinear forecasting for the classification of natural time series. Philosophical Transactions of the Royal Society of London. Series A: Physical and Engineering Sciences, 348(1688), 477-495.
"""
function smap(x::Dataset; θ = 1.0, k = 1, trainees = 1:length(x) ÷ 2, predictees = length(x) ÷ 2 + 1:length(x)-k, metric = Euclidean(), r = 0)
    Mtrain = x[trainees]; nₜ = length(trainees)
    d = size(x, 2)
    nq = length(predictees)

    x̂s = zeros(nq)
    x̂s_truths = zeros(nq)

    # We really only need distances between each query point qᵢ and the points in `Mtrain`, 
    # and we don't need them sorted, so the below approach is not ideal. The direct approach
    # is commented out in the loop below, but that approach doesn't take into consideration
    # the Theiler window (exclude point in a time radius of `r` around each query).
    # Doing a neighbor search with `nₜ` neighbors seems cumbersome, but makes it easier to 
    # use the Theiler window.
    #tree, theiler = KDTree(x[predictees]), Theiler(r)
    #query_neighbors, query_dists = bulksearch(tree, Mtrain, NeighborNumber(nq-k), theiler)
    
    # Make predictions
    for i = 1:nq-k - 1
        # Select the query point
        query = predictees[i]
        xᵢ = x[query]
        
        # Construct weight matrix. These weights take into account the 
        # distance from the query point to all other training points.
        dists = [evaluate(metric, xᵢ, xₜ) for xₜ in Mtrain]
        #dists = query_dists[i]
        
        d̄ = mean(dists)        
        wts = [exp(-θ*d / d̄) for d in dists]
        W = Diagonal(wts[1:nₜ - k])

        # We want to train a model that predicts k time steps into the future
        # For that, each row/point in `A` must be time-shifted `k`` steps 
        # relative to the corresponding entry in `b`. Weights are applied after 
        # selecting the points.
        b = W * Mtrain[1+k:end, 1]
        A = W * [repeat([1], nₜ - k) Matrix(Mtrain[1:nₜ - k])]

        # Least squares solution
        c = A \ b 

        # Eq. 3.1 in Sugihara (1994). Here, we predict x̂ᵢ using the 
        # model obtained from the least squares solution.
        # c[1] = constant, c[2:end] are the coefficients. 
        x̂ᵢ = c[1] + sum(xᵢ .* c[2:end])

        # Locate the corresponding grouth truth
        x̂ᵢ_truth = x[query + k][1]

        x̂s[i] = x̂ᵢ
        x̂s_truths[i] = x̂ᵢ_truth
    end
    
    return x̂s, x̂s_truths
end