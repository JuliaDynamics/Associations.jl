import TimeseriesSurrogates.randomshuffle
import TimeseriesSurrogates.randomphases
import TimeseriesSurrogates.randomamplitudes
import TimeseriesSurrogates.aaft
import TimeseriesSurrogates.iaaft
using DynamicalSystemsBase.Dataset

##############
# Embeddings
##############
"""
    randomshuffle(E::AbstractEmbedding; cols = 1:size(E.points, 1))

Random shuffle surrogate of an embedding. `cols` controls which
variables of the embedding are shuffled.
"""
function randomshuffle(E::AbstractEmbedding; cols = 1:size(E.points, 1))
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
    randomphases(E::AbstractEmbedding; cols = 1:size(E.points, 1))

Random phases Fourier surrogate of an embedding. `cols` controls which
variables of the embedding are shuffled.
"""
function randomphases(E::AbstractEmbedding; cols = 1:size(E.points, 1))
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
    randomamplitudes(E::AbstractEmbedding; cols = 1:size(E.points, 1))

Random amplitude Fourier surrogate of an embedding. `cols` controls which
variables of the embedding are shuffled.
"""
function randomamplitudes(E::AbstractEmbedding; cols = 1:size(E.points, 1))
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
    randomamplitudes(E::AbstractEmbedding; cols = 1:size(E.points, 1))

Amplitude-adjusted Fourier transform (AAFT) surrogate of an embedding. `cols`
controls which variables of the embedding are shuffled.
"""
function aaft(E::AbstractEmbedding; cols = 1:size(E.points, 1))
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
    randomamplitudes(E::AbstractEmbedding; cols = 1:size(E.points, 1))

Iterated amplitude-adjusted Fourier transform (IAAFT) surrogate of an embedding.
`cols` controls which variables of the embedding are shuffled.
"""
function iaaft(E::AbstractEmbedding; cols = 1:size(E.points, 1))
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


##############
# Datasets
##############

"""
    randomshuffle(d::Dataset; cols = 1:size(d, 2))

Random shuffle surrogate of a Dataset.
`cols` controls which variables of the embedding are shuffled.
"""
function randomshuffle(d::Dataset; cols = 1:size(d, 2))
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
    randomphases(d::Dataset; cols = 1:size(d, 2))

Random phases Fourier surrogate of a Dataset.
`cols` controls which variables of the embedding are shuffled.
"""
function randomphases(d::Dataset; cols = 1:size(d, 2))
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
    randomamplitudes(d::Dataset; cols = 1:size(d, 2))

Random amplitude Fourier surrogate of a Dataset.
`cols` controls which variables of the embedding are shuffled.
"""
function randomamplitudes(d::Dataset; cols = 1:size(d, 2))
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
    aaft(d::Dataset; cols = 1:size(d, 2))

Amplitude-adjusted Fourier transform (IAAFT) surrogate of a Dataset.
`cols` controls which variables of the embedding are shuffled.
"""
function aaft(d::Dataset; cols = 1:size(d, 2))
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
    iaaft(d::Dataset; cols = 1:size(d, 2))

Iterated amplitude-adjusted Fourier transform (IAAFT) surrogate of a Dataset.
`cols` controls which variables of the embedding are shuffled.
"""
function iaaft(d::Dataset; cols = 1:size(d, 2))
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




export randomshuffle, randomphases, randomamplitudes, aaft, iaaft
