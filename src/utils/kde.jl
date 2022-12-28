import Entropies: ball_volume
export densities_at_points
export MultivariateKernel, NormalIsotropic, Epanechnikov
export BandwidthRule, Silverman, DiksFang
export bandwidth

abstract type MultivariateKernel end
"""
    Parzen <: MultivariateKernel
    Parzen(h)

The Parzen kernel. For a given `d`-dimensional dataset `x` and a query point `xᵢ`, it
computes the number of points within a hypercube of radius `h` centered on `pᵢ`.

To use it, do `p = Parzen(0.2); p(x, xᵢ)`.
"""
struct Parzen{B} <: MultivariateKernel; h::B; end
#(p::Parzen)(x::AbstractDataset, xᵢ, h) = (p .- xᵢ for p in D) ./ p.h

using SpecialFunctions
using LinearAlgebra

ball_volume(dim::Int, r) = π^(dim / 2) / gamma(dim / 2 +1) * r^dim

struct NormalIsotropic <: MultivariateKernel end
function (k::NormalIsotropic)(x::SVector{D, T}) where {D, T}
    (1/(2π)^(-D/2)) * exp(-0.5 * transpose(x) * x)
end

# -----------------------------------------
# The following kernels are given in
# Statistics 240 Lecture Notes
# P.B. Stark www.stat.berkeley.edu/∼stark/index.html
# https://www.stat.berkeley.edu/~stark/Teach/S240/Notes/ch10.pdf
struct Epanechnikov <: MultivariateKernel end
function (k::Epanechnikov)(x::SVector{D, T}) where {D, T}
    (D + 2) / (2 * ball_volume(D, 1.0)) * (1 - transpose(x)*x)
end

struct Kernel2 <: MultivariateKernel end
function (k::Kernel2)(x::SVector{D, T}) where {D, T}
    (3 / π) * (1 - transpose(x) * x)^2
end

struct Kernel3 <: MultivariateKernel end
function (k::Kernel3)(x::SVector{D, T}) where {D, T}
    (4 / π) * (1 - transpose(x) * x)^3
end
# -----------------------------------------

"""
    densities_at_points(kernel::MultivariateKernel, x::AbstractDataset, bandwidth)

Compute the local densities at each point `xᵢ ∈ x` using the given multivariate `kernel`
and `bandwidth`.
"""
function densities_at_points(kernel::MultivariateKernel, x::AbstractDataset, bandwidth)
    ρs = [density_at_point(kernel, x, bandwidth, xᵢ, i) for (i, xᵢ) in enumerate(x)]
end

function density_at_point(kernel, x::AbstractDataset{D}, bandwidth, xᵢ, i) where D
    ρᵢ = 0.0
    @inbounds for j in eachindex(x)
        if j != i
            ρᵢ += kernel((xᵢ .- x[j]) ./ bandwidth)
        end
    end
    ρᵢ /= ((length(x) - 1)*bandwidth)^D
end



"""
    probability(k::MultivariateKernel, data::AbstractDataset, xᵢ; h) → p̂(xᵢ)

Compute `p̂(xᵢ)`, the kernel density estimate of the probability `p(xᵢ)`, given some
(multivariate) `data` and the query point `xᵢ`.

This is fast if `xᵢ` is an `SVector`.
"""
function probability(kernel::MultivariateKernel, data::AbstractDataset{D, T}, xᵢ;
        h::Real = silvermans_rule(data)) where {D, T}
    n = length(data)
    return 1 / (n * h^D) * sum(kernel((x - p) / h) for p in data)
end

"""
The supertype for all kernel density bandwidth rules.
"""
abstract type BandwidthRule end

"""
    bandwidth(heuristic::BandwidthRule, x::AbstractDataset)

Compute the bandwidth for a kernel density estimator for the input data `x` using
the given `heuristic`.

## Supported heuristics

- [`Silverman`](@ref)
- [`DiksFang`](@ref)
"""
function bandwidth end

"""
    DiksFang <: BandwidthRule
    DickFang(c = 4.8)

A rule of thumb giving the bandwidth ``h`` for kernel density estimation of entropy
as ``h = cn^{-\\dfrac{2}{7}}``, where ``n`` is the number of data points
(Diks & Fang, 2017)[DiksFang2017].

[DiksFang2017]:
    Diks, C., & Fang, H. (2017). Detecting Granger Causality with a Nonparametric
    Information-Based Statistic (No. 17-03). CeNDEF Working Paper.
"""
Base.@kwdef struct DiksFang{C} <: BandwidthRule
    c::C = 4.8
end
function bandwidth(heuristic::DiksFang{C}, x::AbstractDataset) where C
    heuristic.c * length(x)^(-2/7)
end

"""
    Silverman <: BandwidthRule
    Silverman()

A rule of thumb giving the bandwidth ``h`` for kernel density estimation of entropy
following the rules outlined in [this paper](https://hal.archives-ouvertes.fr/hal-00353297/document)

"""
struct Silverman <: BandwidthRule end

function bandwidth(heuristic::Silverman, x::AbstractDataset)
    silvermans_rule(x)
end

function silvermans_rule(data::AbstractDataset{D, T}) where {D, T}
    N = length(data)
    M = vcat(columns(data)...)
    iqr = quantile(M, 0.75) - quantile(M, 0.25)
    return 0.9 * min(std(M), iqr/1.34)*(D*N)^(-1/5)
end
