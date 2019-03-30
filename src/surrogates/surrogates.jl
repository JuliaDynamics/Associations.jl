import TimeseriesSurrogates.randomshuffle
import TimeseriesSurrogates.randomphases
import TimeseriesSurrogates.randomamplitudes
import TimeseriesSurrogates.aaft
import TimeseriesSurrogates.iaaft
import DynamicalSystems.Dataset
import StateSpaceReconstruction: Embeddings


############################
# Multivariate time series
############################

"""
    randomshuffle(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))

Column-wise random shuffle surrogate of an array,
where each column is a scalar-valued time series. `cols` controls which
variables of the embedding are shuffled.
"""
function randomshuffle(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = randomshuffle(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    d_shuffled
end


"""
    randomphases(a::AbstractArray{Number, 2}, cols = 1:size(d, 2))

Column-wise random phases Fourier surrogate of an array,
where each column is a scalar-valued time series. `cols` controls
which variables of the embedding are shuffled.
"""
function randomphases(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = randomphases(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    d_shuffled
end


"""
    randomamplitudes(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))

Column-wise random amplitude Fourier surrogate of an array,
where each column is a scalar-valued time series. `cols` controls
which variables of the embedding are shuffled.
"""
function randomamplitudes(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = randomamplitudes(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    d_shuffled
end


"""
    aaft(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))

Column-wise amplitude-adjusted Fourier transform (AAFT) surrogate of an array,
where each column is a scalar-valued time series.
`cols` controls which variables of the embedding are shuffled.
"""
function aaft(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = aaft(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    d_shuffled
end

"""
    iaaft(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))

Iterated amplitude-adjusted Fourier transform (IAAFT) surrogate of an array,
where each column is a scalar-valued time series.
`cols` controls which variables of the embedding are shuffled.
"""
function iaaft(a::AbstractArray{Number, 2}; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = iaaft(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    Dataset(d_shuffled)
end

##############
# Embeddings
##############
"""
    randomshuffle(E::Embeddings.AbstractEmbedding;
                    cols = 1:size(E.points, 1))

Column-wise random shuffle surrogate of an embedding. `cols` controls which
variables of the embedding are shuffled.
"""
function randomshuffle(E::Embeddings.AbstractEmbedding;
                        cols = 1:size(E.points, 1))
    n_variables = size(E.points, 1)
    E_shuffled = similar(E.points)
    for i = 1:n_variables
        if i ∈ cols
            E_shuffled[i, :] = randomshuffle(E.points[i, :])
        else
            E_shuffled[i, :] = E.points[i, :]
        end
    end

    embed(transpose(E_shuffled))
end

"""
    randomphases(E::Embeddings.AbstractEmbedding; cols = 1:size(E.points, 1))

Column-wise random phases Fourier surrogate of an embedding.
`cols` controls which variables of the embedding are shuffled.
"""
function randomphases(E::Embeddings.AbstractEmbedding;
                        cols = 1:size(E.points, 1))
    n_variables = size(E.points, 1)
    E_shuffled = similar(E.points)
    for i = 1:n_variables
        if i ∈ cols
            E_shuffled[i, :] = randomphases(E.points[i, :])
        else
            E_shuffled[i, :] = E.points[i, :]
        end
    end

    embed(transpose(E_shuffled))
end


"""
    randomamplitudes(E::Embeddings.AbstractEmbedding;
                        cols = 1:size(E.points, 1))

Column-wise random amplitude Fourier surrogate of an embedding. `cols`
controls which variables of the embedding are shuffled.
"""
function randomamplitudes(E::Embeddings.AbstractEmbedding;
                            cols = 1:size(E.points, 1))
    n_variables = size(E.points, 1)
    E_shuffled = similar(E.points)
    for i = 1:n_variables
        if i ∈ cols
            E_shuffled[i, :] = randomamplitudes(E.points[i, :])
        else
            E_shuffled[i, :] = E.points[i, :]
        end
    end

    embed(transpose(E_shuffled))
end


"""
    aaft(E::Embeddings.AbstractEmbedding; cols = 1:size(E.points, 1))

Column-wise amplitude-adjusted Fourier transform (AAFT) surrogate of an
embedding. `cols` controls which variables of the embedding are shuffled.
"""
function aaft(E::Embeddings.AbstractEmbedding; cols = 1:size(E.points, 1))
    n_variables = size(E.points, 1)
    E_shuffled = similar(E.points)
    for i = 1:n_variables
        if i ∈ cols
            E_shuffled[i, :] = aaft(E.points[i, :])
        else
            E_shuffled[i, :] = E.points[i, :]
        end
    end

    embed(transpose(E_shuffled))
end

"""
    iaaft(E::Embeddings.AbstractEmbedding; cols = 1:size(E.points, 1))

Column-wise iterated amplitude-adjusted Fourier transform (IAAFT) surrogate
of an embedding. `cols` controls which variables of the embedding are shuffled.
"""
function iaaft(E::Embeddings.AbstractEmbedding; cols = 1:size(E.points, 1),
				n_maxiter = 200, tol = 1e-6, n_windows = 50)
    n_variables = size(E.points, 1)
    E_shuffled = similar(E.points)
    for i = 1:n_variables
        if i ∈ cols
            E_shuffled[i, :] = iaaft(E.points[i, :],
								n_maxiter = n_maxiter,
								tol = tol,
								n_windows = n_windows)
        else
            E_shuffled[i, :] = E.points[i, :]
        end
    end

    embed(transpose(E_shuffled))
end


##############
# Datasets
##############

"""
    randomshuffle(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))

Column-wise random shuffle surrogate of a Dataset. `cols` controls which
variables of the embedding are shuffled.
"""
function randomshuffle(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = randomshuffle(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    Dataset(d_shuffled)
end


"""
    randomphases(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))

Column-wise random phases Fourier surrogate of a Dataset. `cols` controls
which variables of the embedding are shuffled.
"""
function randomphases(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = randomphases(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    Dataset(d_shuffled)
end


"""
    randomamplitudes(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))

Column-wise random amplitude Fourier surrogate of a Dataset. `cols` controls
which variables of the embedding are shuffled.
"""
function randomamplitudes(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = randomamplitudes(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    Dataset(d_shuffled)
end


"""
    aaft(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))

Column-wise amplitude-adjusted Fourier transform (AAFT) surrogate of a Dataset.
`cols` controls which variables of the embedding are shuffled.
"""
function aaft(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = aaft(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    Dataset(d_shuffled)
end

"""
    iaaft(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))

Iterated amplitude-adjusted Fourier transform (IAAFT) surrogate of a Dataset.
`cols` controls which variables of the embedding are shuffled.
"""
function iaaft(d::DynamicalSystemsBase.Dataset; cols = 1:size(d, 2))
    n_variables = size(d, 2)
    d_shuffled = zeros(eltype(d), size(d))
    for i = 1:n_variables
        if i ∈ cols
            d_shuffled[:, i] = iaaft(d[:, i])
        else
            d_shuffled[:, i] = d[:, i]
        end
    end

    Dataset(d_shuffled)
end



##########################
# CustomReconstructions
##########################

randomshuffle(d::CustomReconstruction; cols = 1:size(d, 2)) = 
    randomshuffle(d.reconstructed_pts; cols = cols)

randomphases(d::CustomReconstruction; cols = 1:size(d, 2)) = 
    randomphases(d.reconstructed_pts; cols = cols)

randomamplitudes(d::CustomReconstruction; cols = 1:size(d, 2)) = 
    randomamplitudes(d.reconstructed_pts; cols = cols)

aaft(d::CustomReconstruction; cols = 1:size(d, 2)) = 
    aaft(d.reconstructed_pts; cols = cols)
  
iaaft(d::CustomReconstruction; cols = 1:size(d, 2)) = 
    iaaft(d.reconstructed_pts; cols = cols)

export randomshuffle, randomphases, randomamplitudes, aaft, iaaft
