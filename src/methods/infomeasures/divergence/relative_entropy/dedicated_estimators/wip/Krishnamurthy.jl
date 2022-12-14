#using QuadGK
#quadgk(f, -Inf, Inf, rtol=1e-3)

struct Krishnamurthy <: InformationMeasure

end


# internal helper
abstract type KrishnamurthyEstimate end
struct P̂ <: KrishnamurthyEstimate end

function kernel_estimate(::KrishnamurthyEstimate, x::AbstractDataset)

end

# Kernel density estimates of densities at point `p`
p̂(k::Kernel, X::AbstractDataset, p) = (1 / length(X)) * sum(k(Xᵢ - p) for Xᵢ in X)
q̂(k::Kernel, Y::AbstractDataset, p) = (1 / length(Y)) * sum(k(Yᵢ - p) for Yᵢ in Y)
# DS := data splitting; k := kernel with some bandwith meeting assumption 3 in paper.
p̂DS(k::Kernel, X::AbstractDataset, p) = (2 / length(Y)) * sum(k(Xᵢ - p) for Xᵢ in X)
q̂DS(k::Kernel, Y::AbstractDataset, p) = (2 / length(Y)) * sum(k(Yᵢ - p) for Yᵢ in Y)
