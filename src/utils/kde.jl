import Entropies: ball_volume

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

# https://hal.archives-ouvertes.fr/hal-00353297/document
function silvermans_rule(data::AbstractDataset{D, T}) where {D, T}
    N = length(data)
    M = vcat(columns(data)...)
    iqr = quantile(M, 0.75) - quantile(M, 0.25)
    return 0.9 * min(std(M), iqr/1.34)*(D*N)^(-1/5)
end
