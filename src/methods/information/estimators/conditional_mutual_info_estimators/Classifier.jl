using Random
using Flux
import Flux: params
using Distributions
using Statistics

export ClassifierCMI 

"""
    ClassifierCMI <: ConditionalMutualInformationEstimator
    ClassifierCMI()

A [`CMIShannon`](@ref) estimator.
"""
struct ClassifierCMI{M, TB, TT, Tτ, RNG} <: ConditionalMutualInformationEstimator{M}
    definition::M
    B::TB # Number of outer bootstrap iterations
    T::TT # Number of inner iterations
    τ::Tτ # Clipping constant
    rng::RNG
end

function ClassifierCMI(; definition = CMIShannon(),  B = 10, T = 100, τ = 1e-10, 
        rng = Random.default_rng())
    return ClassifierCMI(definition, B, T, τ, rng)
end


# Helper function to split data into train and test sets
function split_test_train(d::StateSpaceSet, ratio=0.5)
    n = length(d)
    idx = randperm(n)
    train_size = floor(Int, ratio * n)
    train_idx = @view idx[1:train_size]
    test_idx = @view idx[(train_size+1):end]
    return StateSpaceSet([d.data[i] for i in train_idx]), 
           StateSpaceSet([d.data[i] for i in test_idx])
end

function association(est::ClassifierCMI, x::StateSpaceSet, y::StateSpaceSet, z::StateSpaceSet)
    (; definition, B, T, τ, rng) = est
    n = length(x)
    I_estimates = zeros(B)

    for b in 1:B
        # Permute the datasets
        perm = randperm(rng, n)
        x_π = StateSpaceSet([x.data[i] for i in perm])
        y_π = StateSpaceSet([y.data[i] for i in perm])
        z_π = StateSpaceSet([z.data[i] for i in perm])

        # Split the permuted dataset
        mid = div(n, 2)
        D_class_joint = (StateSpaceSet(x_π.data[1:mid]), StateSpaceSet(y_π.data[1:mid]), StateSpaceSet(z_π.data[1:mid]))
        D_gen = (StateSpaceSet(x_π.data[(mid+1):end]), StateSpaceSet(y_π.data[(mid+1):end]), StateSpaceSet(z_π.data[(mid+1):end]))

        # Train a generator (here we use a simple shuffling as a placeholder)
        function G(z::StateSpaceSet)
            return StateSpaceSet(D_gen[2].data[shuffle(rng, 1:length(D_gen[2]))])
        end

        # Generate marginal dataset
        y_prime = G(D_class_joint[3])
        D_class_marg = (D_class_joint[1], y_prime, D_class_joint[3])
        
        I_estimates[b] = classifier_Dkl(D_class_joint, D_class_marg, T, τ)
    end

    return mean(I_estimates)
end

function classifier_Dkl(D_p, D_q, T, τ)
    D_KL_estimates = zeros(T)

    for t in 1:T
        x_train_p, x_eval_p = split_test_train(D_p[1])
        y_train_p, y_eval_p = split_test_train(D_p[2])
        z_train_p, z_eval_p = split_test_train(D_p[3])
        x_train_q, x_eval_q = split_test_train(D_q[1])
        y_train_q, y_eval_q = split_test_train(D_q[2])
        z_train_q, z_eval_q = split_test_train(D_q[3])

        # Combine datasets for training
        X_train_p = hcat(Matrix(x_train_p)', Matrix(y_train_p)', Matrix(z_train_p)')
        X_train_q = hcat(Matrix(x_train_q)', Matrix(y_train_q)', Matrix(z_train_q)')
        X_train = hcat(X_train_p, X_train_q)
        y_train = vcat(ones(size(X_train_p, 2)), zeros(size(X_train_q, 2)))

        # Define and train the classifier
        input_dim = size(X_train, 1)
        model = Chain(
            Dense(input_dim, 64, relu),
            Dense(64, 32, relu),
            Dense(32, 1, σ)
        )

        loss(x, y) = Flux.binarycrossentropy(model(x), y)
        opt = ADAM()
        data = Flux.Data.DataLoader((X_train, y_train'), batchsize=32, shuffle=true)

        # Manual training loop
        epochs = 10
        for epoch in 1:epochs
            Flux.train!(loss, params(model), data, opt)
        end

        # Get predictions
        X_eval_p = hcat(Matrix(x_eval_p)', Matrix(y_eval_p)', Matrix(z_eval_p)')
        X_eval_q = hcat(Matrix(x_eval_q)', Matrix(y_eval_q)', Matrix(z_eval_q)')
        pred_p = clamp.(vec(model(X_eval_p)), τ, 1-τ)
        pred_q = clamp.(vec(model(X_eval_q)), τ, 1-τ)

        # Compute D_KL estimate
        D_KL_estimates[t] = mean(log.(pred_p ./ (1 .- pred_p))) -
                            log(mean(pred_q ./ (1 .- pred_q)))
    end

    return mean(D_KL_estimates)
end

export generate_gaussian_cmi_data
function generate_gaussian_cmi_data(n_samples::Int; rng::AbstractRNG = Random.default_rng())
    # Define the covariance matrix
    ρ_xy = 0.7  # Correlation between X and Y
    ρ_xz = 0.5  # Correlation between X and Z
    ρ_yz = 0.3  # Correlation between Y and Z
    
    Σ = [1.0    ρ_xy   ρ_xz;
         ρ_xy   1.0    ρ_yz;
         ρ_xz   ρ_yz   1.0]
    
    # Create the multivariate normal distribution
    d = MvNormal(zeros(3), Σ)
    
    # Generate samples
    samples = rand(rng, d, n_samples)
    
    # Convert to StateSpaceSets
    x = StateSpaceSet([SVector{1}(samples[1,i]) for i in 1:n_samples])
    y = StateSpaceSet([SVector{1}(samples[2,i]) for i in 1:n_samples])
    z = StateSpaceSet([SVector{1}(samples[3,i]) for i in 1:n_samples])
    
    return x, y, z, Σ
end
export analytical_gaussian_cmi
function analytical_gaussian_cmi(Σ::Matrix{Float64})
    # Extract the relevant parts of the covariance matrix
    σ_x = Σ[1,1]
    σ_y = Σ[2,2]
    σ_z = Σ[3,3]
    σ_xy = Σ[1,2]
    σ_xz = Σ[1,3]
    σ_yz = Σ[2,3]
    
    # Compute the conditional covariance
    Σ_xy_given_z = [σ_x σ_xy; σ_xy σ_y] - [σ_xz; σ_yz] * [σ_xz σ_yz] / σ_z
    
    # Compute the CMI
    cmi = 0.5 * log(det(Σ_xy_given_z[1,1]) * det(Σ_xy_given_z[2,2]) / det(Σ_xy_given_z))
    
    return cmi
end

# # Helper function to split data into train and test sets
# function split_test_train(x, y, z, ratio=0.5)
#     n = size(x, 2)
#     idx = shuffle(1:n)
#     train_idx = idx[1:floor(Int, ratio*n)]
#     test_idx = idx[ceil(Int, ratio*n):end]
#     return (x[:, train_idx], y[:, train_idx], z[:, train_idx]), 
#            (x[:, test_idx], y[:, test_idx], z[:, test_idx])
# end

# function association(est::ClassifierCMI, x, y, z) 
#     (; definition, B, T, τ, rng) = est
#     n = size(x, 2)
#     I_estimates = zeros(B)
    
#     for b in 1:B
#         # Permute the datasets
#         perm = shuffle(rng, 1:n)
#         x_π, y_π, z_π = x[:, perm], y[:, perm], z[:, perm]
        
#         # Split the permuted dataset
#         mid = div(n, 2)
#         D_class_joint = (x_π[:, 1:mid], y_π[:, 1:mid], z_π[:, 1:mid])
#         D_gen = (x_π[:, mid+1:end], y_π[:, mid+1:end], z_π[:, mid+1:end])
        
#         # Train a generator (here we use a simple shuffling as a placeholder)
#         function G(z)
#             return D_gen[2][:, shuffle(rng, 1:size(D_gen[2], 2))]  # Shuffle y values
#         end

#         # Generate marginal dataset
#         y_prime = G(D_class_joint[3])
#         D_class_marg = (D_class_joint[1], y_prime, D_class_joint[3])
#         I_estimates[b] = classifier_Dkl(D_class_joint, D_class_marg, T, τ)
#     end
    
#     return mean(I_estimates)
# end

# function classifier_Dkl(D_p, D_q, T, τ)
#     D_KL_estimates = zeros(T)
    
#     for t in 1:T
#         (x_train_p, y_train_p, z_train_p), (x_eval_p, y_eval_p, z_eval_p) = split_test_train(D_p...)
#         (x_train_q, y_train_q, z_train_q), (x_eval_q, y_eval_q, z_eval_q) = split_test_train(D_q...)
        
#         # Combine datasets for training
#         X_train = hcat(vcat(x_train_p, y_train_p, z_train_p), vcat(x_train_q, y_train_q, z_train_q))
#         y_train = vcat(ones(size(x_train_p, 2)), zeros(size(x_train_q, 2)))
        
#         # Define and train the classifier
#         input_dim = size(X_train, 1)
#         model = Chain(
#             Dense(input_dim, 64, relu),
#             Dense(64, 32, relu),
#             Dense(32, 1, σ)
#         )
        
#         loss(x, y) = Flux.binarycrossentropy(model(x), y)
#         opt = ADAM()
#         data = Flux.Data.DataLoader((X_train, y_train'), batchsize=32, shuffle=true)
        
#         # Manual training loop
#         epochs = 10
#         for epoch in 1:epochs
#             Flux.train!(loss, params(model), data, opt)
#         end
        
#         # Get predictions
#         pred_p = clamp.(vec(model(vcat(x_eval_p, y_eval_p, z_eval_p))), τ, 1-τ)
#         pred_q = clamp.(vec(model(vcat(x_eval_q, y_eval_q, z_eval_q))), τ, 1-τ)
        
#         # Compute D_KL estimate
#         D_KL_estimates[t] = mean(log.(pred_p ./ (1 .- pred_p))) - 
#                             log(mean(pred_q ./ (1 .- pred_q)))
#     end
    
#     return mean(D_KL_estimates)
# end
