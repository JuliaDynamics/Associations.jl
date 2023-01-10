
import ComplexityMeasures: TransferOperator, invariantmeasure, InvariantMeasure, Probabilities
using ComplexityMeasures.GroupSlices
export TransferOperator

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


function _marginal_encodings(encoder::RectangularBinEncoding, x::VectorOrDataset...)
    X = Dataset(Dataset.(x)...)
    bins = [vec(encode_as_tuple(encoder, pt))' for pt in X]
    joint_bins = reduce(vcat, bins)
    idxs = size.(x, 2) #each input can have different dimensions
    s = 1
    encodings = Vector{Dataset}(undef, length(idxs))
    for (i, cidx) in enumerate(idxs)
        variable_subset = s:(s + cidx - 1)
        s += cidx
        y = @views joint_bins[:, variable_subset]
        encodings[i] = Dataset(y)
    end

    return encodings
end

function transferentropy(
        measure::TransferEntropy,
        est::TransferOperator{<:RectangularBinning}, x...)
    e = measure.e
    joint_pts, vars, τs, js = te_embed(measure.embedding, x...)
    iv = invariantmeasure(joint_pts, est.binning)

    # TODO: this needs to be done more cleverly in ComplexityMeasures.jl, so we don't
    # need to do the conversion twice. We should explicitly store the bin indices for all
    # marginals, not a single encoding integer for each bin. Otherwise, we can't
    # properly subset marginals here and relate them to the approximated invariant measure.
    # The bins visited by the orbit are
    visited_bins_coordinates = Dataset(decode.(Ref(iv.to.encoder), iv.to.bins))
    unique_visited_bins = _marginal_encodings(iv.to.encoder, visited_bins_coordinates)[1]

    # # The subset of visited bins with nonzero measure
    inds_non0measure = findall(iv.ρ .> 0)
    positive_measure_bins = unique_visited_bins[inds_non0measure]

    # Estimate marginal probability distributions from joint measure
    cols_ST = [vars.S; vars.T; vars.C]
    cols_TTf = [vars.Tf; vars.T; vars.C]
    cols_T = [vars.T; vars.C]
    p_T  = marginal_probs_from_μ(cols_T, positive_measure_bins, iv, inds_non0measure)
    p_ST = marginal_probs_from_μ(cols_ST, positive_measure_bins, iv, inds_non0measure)
    p_TTf = marginal_probs_from_μ(cols_TTf, positive_measure_bins, iv, inds_non0measure)
    p_joint = iv.ρ[inds_non0measure]

    te = entropy(e, Probabilities(p_ST)) +
        entropy(e, Probabilities(p_TTf)) -
        entropy(e, Probabilities(p_T)) -
        entropy(e, Probabilities(p_joint))
end

transferentropy(est::TransferOperator{<:RectangularBinning}, s, t; kwargs...) =
    transferentropy(Shannon(; base), est, s, t; kwargs...)
transferentropy(est::TransferOperator{<:RectangularBinning}, s, t, c; kwargs...) =
    transferentropy(Shannon(; base), est, s, t, c; kwargs...)
