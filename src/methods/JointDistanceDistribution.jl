using HypothesisTests
using Distances: Euclidean
using Distances: pairwise
using DelayEmbeddings: genembed

import HypothesisTests: OneSampleTTest, pvalue
export OneSampleTTest, pvalue
export JointDistanceDistribution
export jdd

function normalise_minmax(x::T, vmin, vmax) where T
    if x == zero(T)
        return zero(T)
    else
        return (x - vmin)/(vmax - vmin)
    end
end

"""
    JointDistanceDistribution <: AssociationMeasure end
    JointDistanceDistribution(; metric = Euclidean(), B = 10, D = 2, τ = 1, μ = 0.0)

The joint distance distribution (JDD) measure (Amigó & Hirata, 2018)[^Amigo2018].

## Keyword arguments

- **`distance_metric::Metric`**: An instance of a valid distance metric from `Distances.jl`.
    Defaults to `Euclidean()`.
- **`B::Int`**: The number of equidistant subintervals to divide the interval `[0, 1]` into
    when comparing the normalised distances.
- **`D::Int`**: Embedding dimension.
- **`τ::Int`**: Embedding delay.
- **`μ`**: The hypothetical mean value of the joint distance distribution if there
    is no coupling between `x` and `y` (default is `μ = 0.0`).

## Independence testing

If used with [`independence`](@ref), then a `HypothesisTests.OneSampleTTest`
is returned. Use `pvalue` with keyword `tail = :right` on the result to get the
p-value at 95% confidence.

[^Amigo2018]:
    Amigó, José M., and Yoshito Hirata. "Detecting directional couplings from multivariate
    flows by the joint distance distribution." Chaos: An Interdisciplinary Journal of
    Nonlinear Science 28.7 (2018): 075302.
"""
Base.@kwdef struct JointDistanceDistribution{M, T} <: AssociationMeasure
    metric::M = Euclidean()
    B::Int = 10
    D::Int = 2
    τ::Int = 1
    μ::T = 0.0
end

# The convenience wrapper `jdd`` is in deprecations folder for now.

function estimate(measure::JointDistanceDistribution, est::Nothing, source, target)
    return estimate(measure, source, target)
end

# # Internal method for compatibility with independence tests.
# function estimate(measure::JointDistanceDistribution, source, target)
#     (; metric, B, D, τ) = measure
#     js = ([1 for i = 1:D]...,)
#     τs = (collect(0:-τ:-(D-1)*τ)...,)
#     Ex = genembed(source, τs, js)
#     Ey = genembed(target, τs, js)
#     Mx = Matrix(Ex)
#     My = Matrix(Ey)

#     npts = length(Ex)
#     Dx = pairwise(metric, Mx, Mx, dims = 1)
#     Dy = pairwise(metric, My, My, dims = 1)

#     # Normalise the distances to the interval [0, 1]
#     Dx_min = minimum(Dx[Dx .> 0])
#     Dy_min = minimum(Dy[Dy .> 0])
#     Dx_max = maximum(Dx[Dx .> 0])
#     Dy_max = maximum(Dy[Dy .> 0])

#     Dx_norm = zeros(Float64, size(Dx))
#     Dy_norm = zeros(Float64, size(Dy))

#     for i in LinearIndices(Dx[Dx .> 0])
#         Dx_norm[i] = normalise_minmax(Dx[i], Dx_min, Dx_max)
#     end

#     for i in LinearIndices(Dy[Dy .> 0])
#         Dy_norm[i] = normalise_minmax(Dy[i], Dy_min, Dy_max)
#     end

#     mins_δ_yi_yj = fill(2.0, 2*B)

#     for (k, b) in enumerate(1:2*B)
#         bmin = (b-1)/(2*B)
#         bmax = b/(2*B)

#         # Find the indices of all pairs (i, j) in Dx whose distances fall inside the interval Ib
#         #idxs_Dxs_in_Ib = findall(Dx[])

#         # We don't need to store any indices or distances explicitly, but only
#         # keep track of whether a smaller distance has has been detected.
#         # The maximum possible distance after normalisation is 1.0, so this
#         # value can only decrease as we update.
#         min_δ_yi_yj = 1.0

#         for i = 1:npts
#             for j = (i+1):npts
#                 δ_xi_xj = Dx_norm[i, j]
#                 if bmin < δ_xi_xj <= bmax
#                     δ_yi_yj = Dy_norm[i, j]
#                     if δ_yi_yj < min_δ_yi_yj
#                         min_δ_yi_yj = δ_yi_yj
#                     end
#                 end
#             end
#         end
#         mins_δ_yi_yj[k] = min_δ_yi_yj
#     end

#     Δjdd = [mins_δ_yi_yj[B + i] - mins_δ_yi_yj[i] for i in 1:B]

#     return Δjdd
# end


# Internal method for compatibility with independence tests.
function estimate(measure::JointDistanceDistribution, source, target)
    (; metric, B, D, τ) = measure
    js = ([1 for i = 1:D]...,)
    τs = (collect(0:-τ:-(D-1)*τ)...,)
    Ex = genembed(source, τs, js)
    Ey = genembed(target, τs, js)
    # Mx = Matrix(Ex)
    # My = Matrix(Ey)

    npts = length(Ex)
    Dx = pairwise(metric, Ex.data)
    Dy = pairwise(metric, Ey.data)

    # Normalise the distances to the interval [0, 1]
    fDx = filter(dᵢ -> dᵢ .> 0, Dx)
    fDy = filter(dᵢ -> dᵢ .> 0, Dy)
    Dx_min, Dx_max = minimum(fDx), maximum(fDx)
    Dy_min, Dy_max = minimum(fDy), maximum(fDy)

    # # TODO: the runtime of these loops can be reduced to half by only considering the
    # # lower triangular part of the distance matrix and mirroring it.
    for i in eachindex(Dx)
        Dx[i] = normalise_minmax(Dx[i], Dx_min, Dx_max)
    end

    for i in eachindex(Dy)
        Dy[i] = normalise_minmax(Dy[i], Dy_min, Dy_max)
    end



    δs = fill(2.0, 2 * B)

    for (k, b) in enumerate(1:(2 * B))
        bmin = (b - 1) / 2B
        bmax = b / 2B

        δ̃min = jdd_step3(Dx, Dy, bmin, bmax)
        δs[k] = δ̃min
    end

    Δjdd = [δs[B + i] - δs[i] for i in 1:B]

    #return Δjdd

    #     # Find the indices of all pairs (i, j) in Dx whose distances fall inside the interval Ib
    #     #idxs_Dxs_in_Ib = findall(Dx[])

    #     # We don't need to store any indices or distances explicitly, but only
    #     # keep track of whether a smaller distance has has been detected.
    #     # The maximum possible distance after normalisation is 1.0, so this
    #     # value can only decrease as we update.
    #     min_δ_yi_yj = 1.0

    #     for i = 1:npts
    #         for j = (i+1):npts
    #             δ_xi_xj = Dx[i, j]
    #             if bmin < δ_xi_xj <= bmax
    #                 δ_yi_yj = Dy[i, j]
    #                 if δ_yi_yj < min_δ_yi_yj
    #                     min_δ_yi_yj = δ_yi_yj
    #                 end
    #             end
    #         end
    #     end
    #     mins_δ_yi_yj[k] = min_δ_yi_yj
    # end

end

function jdd_step3(Dx, Dy, bmin, bmax)
    Dy_min = 1.0 # by definition, 1.0 is the largest normalized distance.

    for (k, Dxₖ) in enumerate(Dx)
        if bmin <= Dxₖ < bmax
            Dyₖ = Dy[k]
            if Dyₖ < Dy_min
                Dy_min = Dyₖ
            end
        end
    end
    return Dy_min
end
