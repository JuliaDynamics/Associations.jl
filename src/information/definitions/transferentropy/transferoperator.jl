
import ComplexityMeasures: TransferOperator, invariantmeasure, InvariantMeasure, Probabilities
using ComplexityMeasures.GroupSlices
export TransferOperator

using ComplexityMeasures: Probabilities

"""
	marginal_indices(x)

Returns a column vector `v` with the same number of elements as there are unique
elements in `x`. `v[i]` is the indices of elements in `x` matching `v[i]`.

For example, if the third unique element in `x`, and the element `u₃ = unique(x)[3]`
appears four times in `x`, then `v[3]` is a vector of four integers indicating the
position of the elements matching `u₃`.
"""
function marginal_indices(visited_bins, selected_axes)
    marginal_pts = [x[selected_axes] for x in visited_bins]
    groupinds(groupslices(marginal_pts))
end

"""
    marginal_probs_from_μ(seleced_axes, visited_bins, iv::InvariantMeasure, inds_non0measure)

Estimate marginal probabilities from a pre-computed invariant measure, given a set
of visited bins, an invariant measure and the indices of the positive-measure bins.
The indices in `selected_axes` determines which marginals are selected.
"""
function marginal_probs_from_μ(seleced_axes, visited_bins, iv::InvariantMeasure, inds_non0measure)

    marginal_inds::Vector{Vector{Int}} =
        marginal_indices(visited_bins, seleced_axes)

    # When the invariant measure over the joint space is already known, we don't
    # need to estimate histograms. We simply sum over the nonzero entries of the
    # (already estimated) invariant distribution `iv` in the marginal space
    # (whose indices are given by `seleced_axes`).
    μpos = iv.ρ[inds_non0measure]
    marginal = zeros(Float64, length(marginal_inds))
    @inbounds for i in eachindex(marginal_inds)
        marginal[i] = sum(μpos[marginal_inds[i]])
    end
    return marginal
end


function _marginal_encodings(encoder::RectangularBinEncoding, x::VectorOrStateSpaceSet...)
    X = StateSpaceSet(StateSpaceSet.(x)...)
    bins = [vec(encode_as_tuple(encoder, pt))' for pt in X]
    joint_bins = reduce(vcat, bins)
    idxs = size.(x, 2) #each input can have different dimensions
    s = 1
    encodings = Vector{StateSpaceSet}(undef, length(idxs))
    for (i, cidx) in enumerate(idxs)
        variable_subset = s:(s + cidx - 1)
        s += cidx
        y = @views joint_bins[:, variable_subset]
        encodings[i] = StateSpaceSet(y)
    end

    return encodings
end

# Only works for `RelativeAmount`, because probabilities are obtained from the 
# transfer operator.
function h4_marginal_probs(
        est::EntropyDecomposition{
            <:TransferEntropy, 
            <:DiscreteInfoEstimator, 
            <:TransferOperator, 
            <:RelativeAmount
        },
        x...
    )
    if !est.discretization.binning.precise
        throw(ArgumentError("Please supply a binning with `precise == true`, otherwise points may end up outside the binning."))
    end
    joint_pts, vars, τs, js = te_embed(est.definition.embedding, x...)
    iv = invariantmeasure(joint_pts, est.discretization.binning)

    # TODO: this needs to be done more cleverly in ComplexityMeasures.jl, so we don't
    # need to do the conversion twice. We should explicitly store the bin indices for all
    # marginals, not a single encoding integer for each bin. Otherwise, we can't
    # properly subset marginals here and relate them to the approximated invariant measure.
    encoding = iv.to.encoding
    visited_bins_coordinates = StateSpaceSet(decode.(Ref(encoding), iv.to.bins))
    unique_visited_bins = _marginal_encodings(iv.to.encoding, visited_bins_coordinates)[1]

    # # The subset of visited bins with nonzero measure
    inds_non0measure = findall(iv.ρ .> 0)
    positive_measure_bins = unique_visited_bins[inds_non0measure]

    # Estimate marginal probability distributions from joint measure
    cols_STC = [vars.S; vars.T; vars.C]
    cols_T⁺TC = [vars.Tf; vars.T; vars.C]
    cols_TC = [vars.T; vars.C]
    pTC  = marginal_probs_from_μ(cols_TC, positive_measure_bins, iv, inds_non0measure)
    pSTC = marginal_probs_from_μ(cols_STC, positive_measure_bins, iv, inds_non0measure)
    pT⁺TC = marginal_probs_from_μ(cols_T⁺TC, positive_measure_bins, iv, inds_non0measure)
    pST⁺TC = iv.ρ[inds_non0measure]

    return Probabilities(pTC), 
        Probabilities(pSTC), 
        Probabilities(pT⁺TC), 
        Probabilities(pST⁺TC)
end

function information(
        est::EntropyDecomposition{
            <:TransferEntropy, 
            <:DiscreteInfoEstimator, 
            <:TransferOperator, 
            <:RelativeAmount
        },
        x...)
    # If a conditional input (x[3]) is not provided, then C is just a 0-dimensional
    # StateSpaceSet. The horizontal concatenation of C with T then just returns T.
    # We therefore don't need separate methods for the conditional and non-conditional
    # cases.
    pTC, pSTC, pT⁺TC, pST⁺TC = h4_marginal_probs(est, x...)
    cmi_est = convert_to_cmi_estimator(est)
    h_est = estimator_with_overridden_parameters(cmi_est.definition, cmi_est.est)

    # Estimate by letting TE(s -> t | c) := I(t⁺; s⁻ | t⁻, c⁻).
    hSTC =  information(h_est, pSTC)
    hT⁺TC = information(h_est, pT⁺TC)
    hTC = information(h_est, pTC)
    hST⁺TC = information(h_est, pST⁺TC)
    te = hT⁺TC - hTC - hST⁺TC + hSTC
    return te 
end