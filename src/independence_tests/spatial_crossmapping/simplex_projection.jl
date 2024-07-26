using DelayEmbeddings, Neighborhood

export simplex_predictions

"""
    simplex_predictions(x::AbstractVector, k; τ = 1, d = 2, p = 1, 
        training = 1:length(x) ÷ 2, prediction = length(x) ÷ 2 + 1:length(x)
    ) → Vector{Float64}, Vector{Float64}

From an embedding `Mₓ` computed from the time series `x`, first compute the 
the `k`-step forward simplex projections for each predictee in `Mₓ[prediction]`.
This is done by locating, for each predictee, its `d+1` nearest neighbors in 
`Mₓ[training[1:end-p]]`, and projecting each neighbor `p` steps forward in time. 
Then, for each scalar value `x[i]` where `i ∈ prediction`, compute the prediction 
`x̃[i]` as an exponentially weighted average of the forward-projected neighbors 
(there are `d+1` neighbors, so these form a simplex around the predictee).
If no third argument is given, then two vectors are returned: the observed values and 
the predicted values. If a `correspondence_measure` is given as the third argument, 
then return the correspondence between observed and predicted values.

## Example
```julia
L = 2000
xs = 0.0:0.05:L
ts = sin.(xs) .+ rand(xs) ./ 0.2
# Compute the 3-step forward in time predictions and compute the correlation between
# observed and predicted values.
simplex_predictions(x, 1, Statistics.cor)
```

## Implementation details

This is an implementation of the simplex projection method from Sugihara and May (1990)[^Sugihara1990]. The original 
paper doesn not provide sufficiently detailed pseudocode for implementation, so the algorithm here 
is based on Ye et al. (2015)[^Ye2015].

[^Sugihara1990]: Sugihara, George, and Robert M. May. "Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series." Nature 344.6268 (1990): 734-741.
[^Ye2015]: Ye, H., Beamish, R. J., Glaser, S. M., Grant, S. C., Hsieh, C. H., Richards, L. J., ... & Sugihara, G. (2015). Equation-free mechanistic ecosystem forecasting using empirical dynamic modeling. Proceedings of the National Academy of Sciences, 112(13), E1569-E1576.
"""
function simplex_predictions(x, k; τ = 1, d = 2, 
        training = 1:length(x) ÷ 2, 
        # Predictees need somewhere to go when projected forward in time, 
        # so make sure that we exclude predictees that are among the the 
        # last τ*(d-1) - k points (these would result in simplices for 
        # which one of the vertices has no target `k` steps in the future)
        prediction = length(x) ÷ 2 + 1:length(x)-τ*(d-1)-k
    )

    # Generate embedding and separate the embedding points into a training set and 
    # a set of predictees. 
    τs = 0:-τ:-τ*(d-1)
    Mₓ = genembed(x, τs)
    predictees = Mₓ[prediction]

    # TODO: what happens if there is overlap between training and prediction sets?
    #if isempty(training ∩ prediction)
    # For now, the training set excludes the `k` last points, so that 
    # when later finding nearest neighbors, these neighbors have "somewhere to go" 
    # when projected forward in time. This can be done more intelligently, because
    # some library points are lost, and when using overlapping libraries, the 
    # indexing is not that simple. For now, just use non-overlapping libraries.
    theiler = Theiler(0)
    trainees = Mₓ[training[1:end-k]]

    tree = KDTree(trainees)
    predictee_neighbors, predictee_dists = bulksearch(tree, predictees, NeighborNumber(d+1), theiler)

    #These are re-used
    u = zeros(d + 1)
    w = zeros(d + 1)

    # Predict scalar forward projections. We could also predict d+1-dimensional
    # embedding states, but here we do scalars.
    Lₚ = length(predictees)
    x̃ = zeros(Lₚ)

    for i = 1:Lₚ
        dists = predictee_dists[i]
        dist₁ = dists[1]
        @inbounds for j = 1:d+1
            w[j] = exp(-dists[j] / dist₁)
        end
        s = sum(w)

        # Compute the p-step forward prediction of pᵢ.
        # Select from Mₓ directly, because it also contains the points left 
        # out when defining `training`
        projᵢ = x[predictee_neighbors[i] .+ k]

        # The predicted value is the weighted average of the forward-projected
        # states. 
        for j = 1:d+1
            x̃[i] += w[j]*projᵢ[j] / s
        end
    end

    return x[prediction .+ k], x̃    
end