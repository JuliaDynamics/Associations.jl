using DelayEmbeddings

export delay_simplex

"""
    delay_simplex(x, τ; ds = 2:10, ks = 1:10) → ρs

Determine the optimal embedding dimension for `x` based on the simplex projection algorithm 
from Sugihara & May (1990)[^Sugihara1990]. 

For each `d ∈ ds`, we compute the correlation between observed and predicted values 
for different prediction times `ks`, and average the correlation coefficients. The 
embedding dimension for which the average correlation is highest is taken as the optimal 
dimension. The embedding delay `τ` is given as a positive number.

Returns the prediction skills `ρs` - one `ρ` for each `d ∈ ds`.

Note: the library/training and prediction sets are automatically taken as the first and 
second halves of the data, respectively. This convenience method does not allow tuning 
the libraries further.

[^Sugihara1990]: Sugihara, George, and Robert M. May. "Nonlinear forecasting as a way of distinguishing chaos from measurement error in time series." Nature 344.6268 (1990): 734-741.
"""
function delay_simplex(x, τ; ds = 2:10, ks = 1:8)
    ρs = zeros(length(ds))
    for (i, d) in enumerate(ds)
        embedding = genembed(x, -τ)

        ρs_k = zeros(length(ks))
        for (j, k) in enumerate(ks)
            train = 1:length(x) ÷ 2 - k
            pred = length(x) ÷ 2 + 1 : length(x) - d*τ

            x̄, x̄_actual = simplex_predictions(x, k; τ = τ, d = d, 
                training = train, 
                prediction = pred)
            ρs_k[j] = cor(x̄, x̄_actual)
        end
        ρs[i] = mean(ρs_k)
    end

    return ρs
end