# FOCI Implementation Pseudocode

using Statistics
# Required imports for independence testing
# using LocalPermutationTest, AzadkiaChatterjeeCoefficient, independence

export FOCI
export FOCISelectedFeatures

"""
    FOCI <: GraphAlgorithm
    FOCI(; unconditional_measure = ChatterjeeCorrelation(),
           conditional_measure = AzadkiaChatterjeeCoefficient(),
           unconditional_test = nothing,
           use_significance_test::Bool = false,
           standardize::Bool = true,
           max_vars_per_step::Int = 1,
           forward_backward::Bool = false,
           α = 0.05)

Feature Ordering by Conditional Independence (FOCI) algorithm for feature selection.

## Parameters
- `unconditional_measure`: Measure for unconditional dependence (typically ChatterjeeCorrelation)
- `conditional_measure`: Measure for conditional dependence (typically AzadkiaChatterjeeCoefficient)  
- `unconditional_test`: Test for unconditional significance (e.g., SurrogateAssociationTest)
- `use_significance_test`: If true, use significance testing; if false, use T ≤ 0 criterion
- `standardize`: Whether to standardize features
- `forward_backward`: Whether to use backward elimination
- `α`: Significance level for tests
"""
Base.@kwdef struct FOCI{U,C,T} <: GraphAlgorithm
    unconditional_measure::U = ChatterjeeCorrelation()
    conditional_measure::C = AzadkiaChatterjeeCoefficient()
    unconditional_test::T = nothing
    use_significance_test::Bool = false
    standardize::Bool = true
    max_vars_per_step::Int = 1
    forward_backward::Bool = false
    α = 0.05
end

"""
    FOCISelectedFeatures

Storage for features selected by FOCI algorithm.
"""
Base.@kwdef mutable struct FOCISelectedFeatures{F,FI}
    all_feature_idxs::Vector{Int}
    selected_features::F = Vector{Vector{Float64}}(undef, 0)
    selected_feature_idxs::FI = Vector{Int}(undef, 0)
end

function infer_graph(alg::FOCI, x; verbose=true)
    print_status(FOCIInfoMessage(); verbose)
    return select_features(alg, x; verbose)
end

function infer_graph(alg::FOCI, x::AbstractStateSpaceSet; verbose=true)
    return infer_graph(alg, columns(x); verbose)
end

"""
    select_features(alg::FOCI, x)

The feature selection step of the [`FOCI`](@ref) algorithm, which identifies the
relevant features for each `xᵢ ∈ x`, assuming that `x` must be integer-indexable, i.e.
`x[i]` yields the `i`-th variable.
"""
function select_features(alg::FOCI, x; verbose=false)
    # Find the relevant features for each variable as response.
    features = [select_features(alg, x, k; verbose) for k in eachindex(x)]
    return features
end

function select_features(alg::FOCI, x, i::Int; verbose=false)
    # Extract response variable and predictors
    y = x[i]  # variable i is the response
    X = [x[j] for j in eachindex(x) if j != i]  # all other variables are predictors

    verbose && printstyled("\nSelecting features for response variable "; color=:default)
    verbose && printstyled("x$i\n"; color=:green)

    # Verify input dimensions
    n = length(y)
    for (j, xⱼ) in enumerate(X)
        if length(xⱼ) != n
            throw(ArgumentError("All predictor variables must have the same length as the response variable. " *
                                "Got length(y) = $n, but length(X[$j]) = $(length(xⱼ))."))
        end
    end

    # Standardize if requested
    if alg.standardize
        X = standardize_features(X)
        y = standardize_feature(y)
    end

    verbose && printstyled("Starting FOCI feature selection\n"; bold=true)

    # Initialize feature storage
    features = FOCISelectedFeatures(all_feature_idxs=collect(1:length(X)))

    ###################################################################
    # Forward search
    ###################################################################
    verbose && printstyled("˧ Forward selection...\n"; color=:blue)

    # Step 1: Find j₁ that maximizes Tₙ(Y, Xⱼ)
    significant_feature = select_first_feature!(alg, features, X, y; verbose)

    if significant_feature
        # Step 2: Continue selecting features that maximize conditional dependence
        significant_cond = true
        while significant_cond && length(features.selected_features) < length(X)
            significant_cond = select_next_feature!(alg, features, X, y; verbose)
        end

        ###################################################################
        # Backward elimination (optional)
        ###################################################################
        if alg.forward_backward
            verbose && printstyled("˧ Backward elimination...\n"; color=:blue)
            backward_eliminate!(alg, features, X, y; verbose)
        end
    end

    # Map back to original indices (accounting for removed response variable)
    original_indices = [j for j in eachindex(x) if j != i]

    # Convert selected indices back to original variable indices
    if !isempty(features.selected_feature_idxs)
        features.selected_feature_idxs = [original_indices[j] for j in features.selected_feature_idxs]
    end

    return features
end

# Helper function for when you have explicit X and y (not part of exported API)
function select_features_xy(alg::FOCI, X::Vector{Vector{<:Real}}, y::Vector{<:Real}; verbose=false)
    # Verify input dimensions
    n = length(y)
    for (i, xᵢ) in enumerate(X)
        if length(xᵢ) != n
            throw(ArgumentError("All predictor variables must have the same length as the response variable. " *
                                "Got length(y) = $n, but length(X[$i]) = $(length(xᵢ))."))
        end
    end

    # Standardize if requested
    if alg.standardize
        X = standardize_features(X)
        y = standardize_feature(y)
    end

    verbose && printstyled("Starting FOCI feature selection\n"; bold=true)
    verbose && print_foci_info()

    # Initialize feature storage
    features = FOCISelectedFeatures(all_feature_idxs=1:length(X))

    ###################################################################
    # Forward search
    ###################################################################
    verbose && printstyled("˧ Forward selection...\n"; color=:blue)

    # Step 1: Find j₁ that maximizes Tₙ(Y, Xⱼ)
    significant_feature = select_first_feature!(alg, features, X, y; verbose)

    if significant_feature
        # Step 2: Continue selecting features that maximize conditional dependence
        significant_cond = true
        while significant_cond && length(features.selected_features) < length(X)
            significant_cond = select_next_feature!(alg, features, X, y; verbose)
        end

        ###################################################################
        # Backward elimination (optional)
        ###################################################################
        if alg.forward_backward
            verbose && printstyled("˧ Backward elimination...\n"; color=:blue)
            backward_eliminate!(alg, features, X, y; verbose)
        end
    end

    return features
end

function select_first_feature!(alg::FOCI, features, X, y; verbose=false)
    # Compute Tₙ(Y, Xⱼ) for all features using unconditional measure
    scores = zeros(length(X))
    for (j, xⱼ) in enumerate(X)
        scores[j] = association(alg.unconditional_measure, y, xⱼ)
    end

    # Find feature with maximum score
    best_idx = argmax(scores)
    best_score = scores[best_idx]

    # Apply stopping criterion based on algorithm settings
    if alg.use_significance_test && !isnothing(alg.unconditional_test)
        # Use significance testing
        if best_score <= 0
            verbose && printstyled("  No positive scores found (all T ≤ 0)\n"; color=:yellow)
            return false
        end

        # Test significance of the best feature
        test_result = independence(alg.unconditional_test, y, X[best_idx])

        if test_result.pvalue >= alg.α
            verbose && printstyled("  Best feature not significant (p = $(round(test_result.pvalue, digits=4)) ≥ α)\n"; color=:yellow)
            return false
        end

        verbose && print_feature_selected(best_idx, best_score, "unconditional", test_result.pvalue)
    else
        # Use original FOCI criterion: Tₙ(Y, Xⱼ₁) ≤ 0
        if best_score <= 0
            verbose && printstyled("  No significant features found (T ≤ 0)\n"; color=:yellow)
            return false
        end

        verbose && print_feature_selected(best_idx, best_score, "unconditional")
    end

    # Add the best feature
    push!(features.selected_features, X[best_idx])
    push!(features.selected_feature_idxs, best_idx)

    return true
end

function select_next_feature!(alg::FOCI, features, X, y; verbose=false)
    # Get indices of remaining (unselected) features
    remaining_idxs = setdiff(features.all_feature_idxs, features.selected_feature_idxs)
    isempty(remaining_idxs) && return false

    # Compute conditional dependence measure and test function
    compute_conditional_score, _ = conditional_score_function(alg, features)

    # Compute Tₙ(Y, Xⱼ | S) for all remaining features
    scores = zeros(length(remaining_idxs))
    for (i, j) in enumerate(remaining_idxs)
        scores[i] = compute_conditional_score(y, X[j])
    end

    # Find feature with maximum conditional score
    best_local_idx = argmax(scores)
    best_score = scores[best_local_idx]
    best_global_idx = remaining_idxs[best_local_idx]

    # Check stopping criterion: Tₙ(Y, Xⱼₖ₊₁ | Xⱼ₁, ..., Xⱼₖ) ≤ 0
    if best_score <= 0
        verbose && printstyled("  No more significant conditional associations found (T ≤ 0)\n"; color=:yellow)
        return false
    end

    # Add the best feature
    push!(features.selected_features, X[best_global_idx])
    push!(features.selected_feature_idxs, best_global_idx)

    verbose && print_feature_selected(best_global_idx, best_score, "conditional")
    return true
end

function backward_eliminate!(alg::FOCI, features, X, y; verbose=false)
    length(features.selected_features) < 2 && return features

    verbose && printstyled("  Starting backward elimination...\n"; color=:blue)
    n_initial = length(features.selected_feature_idxs)

    variable_was_eliminated = true
    q = 0
    while variable_was_eliminated && length(features.selected_feature_idxs) >= 2 && q < n_initial
        q += 1
        variable_was_eliminated = elimination_step!(alg, features, X, y; verbose)
    end

    verbose && printstyled("  Backward elimination completed\n"; color=:blue)
    return features
end

function elimination_step!(alg::FOCI, features, X, y; verbose=false)
    isempty(features.selected_features) && return false

    M = length(features.selected_features)
    selected_features = features.selected_features
    selected_idxs = features.selected_feature_idxs

    for k in 1:M
        # Test if feature k is conditionally independent of y given remaining features
        remaining_feature_idxs = setdiff(1:M, k)

        if isempty(remaining_feature_idxs)
            # Only one feature left, can't eliminate further
            continue
        end

        # Create conditioning set from remaining features
        remaining_features = selected_features[remaining_feature_idxs]
        conditioning_set = StateSpaceSet(remaining_features...)

        # Test: y ⫫ selected_features[k] | remaining_features
        # Using LocalPermutationTest with AzadkiaChatterjeeCoefficient
        independence_test = LocalPermutationTest(
            AzadkiaChatterjeeCoefficient(),
            nshuffles=100
        )

        test_result = independence(independence_test, y, selected_features[k], conditioning_set)

        verbose && print_elimination_test(selected_idxs[k], test_result, alg.α)

        # If conditionally independent (p-value >= α), eliminate this feature
        if test_result.pvalue >= alg.α
            verbose && print_feature_eliminated(selected_idxs[k], test_result.pvalue)

            # Remove feature k
            deleteat!(features.selected_features, k)
            deleteat!(features.selected_feature_idxs, k)
            return true
        end
    end

    return false
end

# Similar to your OCE implementation - create conditional scoring function
function conditional_score_function(alg::FOCI, features::FOCISelectedFeatures)
    if isempty(features.selected_features)
        # Pairwise case (shouldn't happen in select_next_feature!, but for completeness)
        compute_score = (y, xⱼ) -> association(alg.unconditional_measure, y, xⱼ)
        test_independence = nothing  # Not used in FOCI
    else
        # Conditional case
        conditioning_set = StateSpaceSet(features.selected_features...)
        compute_score = (y, xⱼ) -> association(alg.conditional_measure, y, xⱼ, conditioning_set)
        test_independence = nothing  # Not used in FOCI
    end
    return compute_score, test_independence
end

function standardize_features(X::Vector{Vector{T}}) where T
    return [standardize_feature(x) for x in X]
end

function standardize_feature(x::Vector{T}) where T
    μ = mean(x)
    σ = std(x)
    return σ > 0 ? (x .- μ) ./ σ : x
end

#########################################################################################
# Pretty printing (following your OCE style)
#########################################################################################

function print_foci_info()
    printstyled("Feature Ordering by Conditional Independence (FOCI)\n"; bold=true)
    printstyled("Notation:\n"; underline=true, color=:default)
    printstyled("  Tₙ(Y, Xⱼ)     := unconditional dependence score\n"; color=:default)
    printstyled("  Tₙ(Y, Xⱼ | S) := conditional dependence score given selected set S\n"; color=:default)
    printstyled("  Y             := response variable\n"; color=:green)
    printstyled("  Xⱼ            := candidate feature j\n"; color=:blue)
    printstyled("  S             := selected feature set\n"; color=:yellow)
    print("\n")
end

struct FOCIInfoMessage end
function print_status(::FOCIInfoMessage; verbose=true)
    if verbose
        print_foci_info()
    end
end

function print_feature_selected(feature_idx::Int, score::Float64, test_type::String, pvalue::Union{Float64,Nothing}=nothing)
    printstyled("  Selected feature "; color=:default)
    printstyled("X$(feature_idx)"; color=:blue)
    if isnothing(pvalue)
        printstyled(" ($test_type score: $(round(score, digits=4)))\n"; color=:default)
    else
        printstyled(" ($test_type score: $(round(score, digits=4)), p = $(round(pvalue, digits=4)))\n"; color=:default)
    end
end

function print_feature_eliminated(feature_idx::Int, pvalue::Float64)
    printstyled("  Eliminated feature "; color=:default)
    printstyled("X$(feature_idx)"; color=:blue)
    printstyled(" (p-value: $(round(pvalue, digits=4)) ≥ α)\n"; color=:default)
end

function print_elimination_test(feature_idx::Int, test_result, α::Float64)
    dep_symbol = test_result.pvalue >= α ? " ⫫ " : " !⫫ "
    action = test_result.pvalue >= α ? "Removing" : "Keeping"

    printstyled("  Y"; color=:green)
    printstyled(dep_symbol; color=:cyan)
    printstyled("X$(feature_idx)"; color=:blue)
    printstyled(" | remaining features → $action "; color=:cyan)
    printstyled("X$(feature_idx)"; color=:blue)
    printstyled(" (p = $(round(test_result.pvalue, digits=4)))\n"; color=:default)
end

# Return type display
function Base.show(io::IO, ::MIME"text/plain", features::FOCISelectedFeatures)
    printstyled("Selected features: "; color=:default)
    if isempty(features.selected_feature_idxs)
        printstyled("∅\n"; color=:yellow)
    else
        printstyled("{"; color=:default)
        for (i, idx) in enumerate(features.selected_feature_idxs)
            printstyled("X$(idx)"; color=:blue)
            i < length(features.selected_feature_idxs) && printstyled(", "; color=:default)
        end
        printstyled("}\n"; color=:default)
    end
end