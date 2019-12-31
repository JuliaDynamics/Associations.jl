"""
Add observational noise to a dataset from a normal distribution with zero mean 
and standard deviation equivalent to some fraction of the empirical standard 
deviation along each variable of the dataset.
"""
function add_observational_noise!(x::AbstractArray{T, 2}, percent_noise) where T
    
    dims = size(x)
    if percent_noise > 0 
        if dims[1] > dims[2]
            nvars = dims[2]
            npts = dims[1]
            
            for i = 1:nvars
                dist = Normal(0, (percent_noise/100) * std(x[:, i]))
                x[:, i] .+= rand(dist, npts)
            end
        else
            nvars = dims[1]
            npts = dims[2]
            
            for i = 1:nvars
                dist = Normal(0, (percent_noise/100) * std(x[i, :]))
                x[i, :] .+= rand(dist, npts)
            end
        end
    end
    
    return x
end

function add_observational_noise!(data::AbstractDataset{dim, T}, percent_noise) where {dim, T}
    D = dim
    L = length(data)
    noisy_data = zeros(T, L, D)
    if percent_noise > 0 
        σs = [Normal(0, (percent_noise/100)*std(col)) for col in columns(data)]
        for i = 1:D
            noisy_data[:, i] = data[:, i] .+ rand(σs[i], L)
        end
        return Dataset(noisy_data)
    end
    return data
end