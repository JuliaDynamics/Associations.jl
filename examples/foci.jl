# Complete FOCI Example: Conservative Approach with Backward Elimination
# =====================================================================

using Random, Statistics, Printf
Random.seed!(1234)

# Generate synthetic data with known relationships
# ===============================================
println("=== Data Generation ===")
n = 500  # number of observations
p = 8    # number of variables (smaller for clearer analysis)

# Create variables with known dependencies
x = [randn(n) for _ in 1:p]

# Define true relationships
x[2] = x[1] .+ 0.5 .* randn(n)                    # x₂ depends on x₁
x[3] = 2.0 .* x[1] .+ x[2] .+ 0.2 .* randn(n)    # x₃ depends on x₁ and x₂  
x[4] = 0.8 .* x[2] .+ 0.3 .* randn(n)             # x₄ depends on x₂
x[5] = -1.5 .* x[1] .+ 0.4 .* randn(n)            # x₅ depends on x₁
# x₆, x₇, x₈ remain independent noise

println("True relationships:")
println("  x₁: independent")
println("  x₂ ← x₁")
println("  x₃ ← x₁, x₂")
println("  x₄ ← x₂")
println("  x₅ ← x₁")
println("  x₆, x₇, x₈: independent")

# Define the true adjacency matrix for comparison
true_adj = zeros(Int, p, p)
true_adj[1, 2] = 1  # x₁ → x₂
true_adj[1, 3] = 1  # x₁ → x₃
true_adj[2, 3] = 1  # x₂ → x₃
true_adj[2, 4] = 1  # x₂ → x₄
true_adj[1, 5] = 1  # x₁ → x₅

println("\nTrue adjacency matrix (rows → columns):")
display(true_adj)

# FOCI Analysis: Multiple Approaches
# ==================================

println("\n=== FOCI Analysis ===")

# 1. Liberal approach (original FOCI)
println("\n1. Liberal FOCI (T ≤ 0 criterion, no backward elimination)")
alg_liberal = FOCI(
    use_significance_test=false,
    forward_backward=false,
    standardize=true
)

results_liberal = infer_graph(alg_liberal, x; verbose=false)

# 2. Conservative approach with significance testing
println("\n2. Conservative FOCI (significance testing, α = 0.01)")
alg_conservative = FOCI(
    unconditional_test=SurrogateAssociationTest(
        ChatterjeeCorrelation();
        nshuffles=200,
        surrogate=RandomShuffle()
    ),
    use_significance_test=true,
    forward_backward=false,
    α=0.01,  # Stricter significance level
    standardize=true
)

results_conservative = infer_graph(alg_conservative, x; verbose=false)

# 3. Full approach: significance testing + backward elimination
println("\n3. Full FOCI (significance testing + backward elimination, α = 0.01)")
alg_full = FOCI(
    unconditional_test=SurrogateAssociationTest(
        ChatterjeeCorrelation();
        nshuffles=200,
        surrogate=RandomShuffle()
    ),
    use_significance_test=true,
    forward_backward=true,  # Enable backward elimination
    α=0.01,
    standardize=true
)

results_full = infer_graph(alg_full, x; verbose=true)

# Analysis Functions
# =================

function results_to_adjacency(results)
    """Convert FOCI results to adjacency matrix"""
    n = length(results)
    adj = zeros(Int, n, n)

    for (j, features) in enumerate(results)  # j is target variable
        for i in features.selected_feature_idxs  # i predicts j
            adj[i, j] = 1
        end
    end

    return adj
end

function analyze_results(results, method_name, true_adj)
    """Analyze FOCI results against ground truth"""
    println("\n--- $method_name Results ---")

    # Show selected features
    for (i, features) in enumerate(results)
        if !isempty(features.selected_feature_idxs)
            predictors = features.selected_feature_idxs
            println("  x$i ← {$(join(["x$j" for j in predictors], ", "))}")
        else
            println("  x$i ← ∅")
        end
    end

    # Convert to adjacency matrix
    pred_adj = results_to_adjacency(results)

    println("\nPredicted adjacency matrix:")
    display(pred_adj)

    # Calculate performance metrics
    n_vars = size(true_adj, 1)

    # True/False Positives/Negatives
    tp = sum((true_adj .== 1) .& (pred_adj .== 1))  # Correctly identified edges
    fp = sum((true_adj .== 0) .& (pred_adj .== 1))  # False edges
    tn = sum((true_adj .== 0) .& (pred_adj .== 0))  # Correctly identified non-edges  
    fn = sum((true_adj .== 1) .& (pred_adj .== 0))  # Missed edges

    # Calculate metrics
    precision = tp / (tp + fp)
    recall = tp / (tp + fn)
    f1 = 2 * (precision * recall) / (precision + recall)
    accuracy = (tp + tn) / (n_vars * n_vars)

    println("\nPerformance Metrics:")
    println("  True Positives:  $tp")
    println("  False Positives: $fp")
    println("  True Negatives:  $tn")
    println("  False Negatives: $fn")
    println("  Precision: $(round(precision, digits=3))")
    println("  Recall:    $(round(recall, digits=3))")
    println("  F1-Score:  $(round(f1, digits=3))")
    println("  Accuracy:  $(round(accuracy, digits=3))")

    return pred_adj, (tp=tp, fp=fp, tn=tn, fn=fn, precision=precision, recall=recall, f1=f1, accuracy=accuracy)
end

# Comparative Analysis
# ===================

println("\n" * "="^60)
println("COMPARATIVE ANALYSIS")
println("="^60)

# Analyze each approach
adj_liberal, metrics_liberal = analyze_results(results_liberal, "Liberal FOCI", true_adj)
adj_conservative, metrics_conservative = analyze_results(results_conservative, "Conservative FOCI", true_adj)
adj_full, metrics_full = analyze_results(results_full, "Full FOCI", true_adj)

# Summary comparison
println("\n" * "="^60)
println("SUMMARY COMPARISON")
println("="^60)

methods = [
    ("Liberal (T ≤ 0)", metrics_liberal),
    ("Conservative (α=0.01)", metrics_conservative),
    ("Full (α=0.01 + BackElim)", metrics_full)
]

println("Method                    | Precision | Recall | F1-Score | Accuracy | TP | FP | TN | FN")
println("-"^80)
for (name, metrics) in methods
    println(@sprintf("%-25s | %7.3f | %6.3f | %8.3f | %8.3f | %2d | %2d | %2d | %2d",
        name, metrics.precision, metrics.recall, metrics.f1, metrics.accuracy,
        metrics.tp, metrics.fp, metrics.tn, metrics.fn))
end

# Detailed Edge Analysis
println("\n" * "="^60)
println("DETAILED EDGE ANALYSIS")
println("="^60)

function analyze_specific_edges(true_adj, pred_adj, method_name)
    println("\n$method_name - Edge-by-Edge Analysis:")
    n = size(true_adj, 1)

    println("Correctly identified edges:")
    for i in 1:n, j in 1:n
        if true_adj[i, j] == 1 && pred_adj[i, j] == 1
            println("  ✅ x$i → x$j")
        end
    end

    println("Missed edges (False Negatives):")
    for i in 1:n, j in 1:n
        if true_adj[i, j] == 1 && pred_adj[i, j] == 0
            println("  ❌ x$i → x$j (missed)")
        end
    end

    println("False edges (False Positives):")
    for i in 1:n, j in 1:n
        if true_adj[i, j] == 0 && pred_adj[i, j] == 1
            println("  ⚠️  x$i → x$j (false positive)")
        end
    end
end

analyze_specific_edges(true_adj, adj_liberal, "Liberal FOCI")
analyze_specific_edges(true_adj, adj_conservative, "Conservative FOCI")
analyze_specific_edges(true_adj, adj_full, "Full FOCI")

# Recommendations
println("\n" * "="^60)
println("RECOMMENDATIONS")
println("="^60)

best_method = argmax([metrics_liberal.f1, metrics_conservative.f1, metrics_full.f1])
method_names = ["Liberal FOCI", "Conservative FOCI", "Full FOCI"]

println("Best performing method: $(method_names[best_method])")
println("\nMethod characteristics:")
println("• Liberal FOCI: Fast, more relationships detected, higher recall")
println("• Conservative FOCI: Fewer false positives, higher precision")
println("• Full FOCI: Most thorough, backward elimination removes spurious relationships")

println("\nUse cases:")
println("• Exploratory analysis: Liberal FOCI")
println("• Feature selection for prediction: Conservative or Full FOCI")
println("• Causal discovery: Full FOCI (but consider dedicated causal methods)")

# Save results
println("\n" * "="^60)
println("RESULTS SUMMARY")
println("="^60)

println("Final adjacency matrices saved:")
println("\nTrue relationships:")
display(true_adj)
println("\nFull FOCI results (recommended):")
display(adj_full)