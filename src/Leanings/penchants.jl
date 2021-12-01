
export mean_observed_penchant, mean_leaning

using Statistics

# We just assume that the time series have been pre-binned or symbolized
# However, `get_states`` might turn out much more complicated if we decide 
# to implement more complicated cause-effect assignments
get_states(x, y) = unique([x; y])

# Getting the penchants might also become more complex for other 
# complicated cause-effect assignments, so separate it out in a function.
function get_penchants(x, y)
    states = get_states(x, y)

    # All possible 2-permutations (with repetitions) of the unique
    # states in `x ∪ y`.
    return Iterators.product(states, states)
end

function penchant_counts(x, y, l = 1)
    @assert length(x) == length(y)
    n = length(x)
    
    # Leave reasonably many for probability computations. `l` can't be too large.
    @assert l < n ÷ 2

    # Number of time the state y_t = 1 AND x_{t-1} = 1 appears in {x, y}
    ns_ec = Vector{Int}(undef, 0)
    ns_c = Vector{Int}(undef, 0)
    ns_e = Vector{Int}(undef, 0)

    # For each possible penchant, count occurrences.
    penchants = get_penchants(x, y)

    @inbounds for penchant in penchants
        sx, sy = penchant[1], penchant[2]
            
        # Number of times the assumed cause has appeared
        n_c = 0

        # Number of times the assumed effect has appeared
        n_e = 0

        # Number of times the cause AND effect has appeared simulaneously
        n_ec = 0

        for t = l+1:n
            # "we required that the cause must precede the effect". This is the 
            # standard l-assignment
            effect = y[t]
            cause = x[t-l]

            if effect == sx && cause == sy
                n_ec += 1
            end
            

            if effect == sx
                n_e += 1
            end

            if cause == sy
                n_c += 1
            end
        end
        push!(ns_ec, n_ec)
        push!(ns_c, n_c)
        push!(ns_e, n_e)
    end

    return ns_c, ns_e, ns_ec
end


function mean_observed_ρ(ρs, κs, ns_ec, L, weighted = false)
    ρ = 0.0
    ct = 0
    @inbounds for i = 1:length(κs)
        if κs[i] != 0.0
            if weighted
                ρ += ρs[i] * ns_ec[i]
            else
                ρ += ρs[i]
            end
            ct += 1
        end
    end

    if weighted
        ρ / L 
    else 
        return ρ / ct
    end
end

"""
    mean_observed_penchant(x, y, l, weighted = false) → ρ̄

Computes the *mean observed penchant* of the `l`-assignment ``\\{C, E\\} = \\{x, y\\}, \\{x_{t-l}, y_t\\}``.

If `weighted == true`, then compute the *weighted mean observed penchant*.

[^McCrackenWeigel2016]: McCracken, J. M., & Weigel, R. S. (2016). Nonparametric causal inference for bivariate time series. Physical Review E, 93(2), 022207.
"""
function mean_observed_penchant(x, y, l = 1; weighted = false)
    n = length(x)
    ns_c, ns_e, ns_ec = penchant_counts(x, y, l)

    L = n - l # "library length", resulting from shifting the time series l steps
    Ps_e = ns_e ./ L
    Ps_c = ns_c ./ L
    κs = ns_ec ./ ns_c
    ρs = κs.* (1 .+ Ps_c ./ (1 .- Ps_c)) .- (Ps_e) ./ (1 .- Ps_c)
    
    # Mean observed penchant. We can determine which penchants are observed 
    # by looking for zeros in `κs`. If `κs[i] == 0`, then the `i`-th penchant
    # is unobserved. Why? When the penchant is unobserved, there are no 
    # instances of the cause and effect occurring together, so 
    # `P(effect | cause) = P(effect, cause) / P(cause) = 0 / P(cause) = 0` 
    # As explained in section D in the paper, including unobserved penchants
    # would be meaningless for causal inference.
    mean_observed_ρ(ρs, κs, ns_ec, L, weighted)
end

"""

    mean_leaning(x, y, l = 1, weighted = true) → ρ̄

Compute the *mean observed leaning* `ρ̄` (McCracken & Weigel, 2016)[^McCrackenWeigel2016] between scalar time series `x` and `y`,
based on the standard `l`-assignment of the `l`-assignment ``\\{C, E\\} = \\{x, y\\}, \\{x_{t-l}, y_t\\}``.

If `ρ̄ > 0`, then the probability that `x` drives `y` is higher than the probability 
that `y` drives `x`. Vice versa, if `ρ̄ > 0`, then the probability that `y` drives `x` is higher than the probability 
that `x` drives `y`.

If `weighted == true`, then each penchant is weighted by the number of times it appears in 
the data, which yields the *weighted mean observed leaning*. This the default behavior, since 
"...[the weighted mean observed leaning] accounts for the frequency of observed cause-effect pairs 
within the data, which is assumed to be a predictor of causal influence" (McCracken & Weigel, 2016).

This implementation *does not* take care of the tolerance domains discussed in the paper. 
Pre-process your time series, using appropriate binning or symbolization schemes, 
so that `x` and `y` both consists of discrete symbols. 

[^McCrackenWeigel2016]: McCracken, J. M., & Weigel, R. S. (2016). Nonparametric causal inference for bivariate time series. Physical Review E, 93(2), 022207.
"""
function mean_leaning(x, y, l = 1, weighted = true)
    return mean_observed_penchant(x, y, l, weighted = weighted) - 
        mean_observed_penchant(y, x, l, weighted = weighted)
end