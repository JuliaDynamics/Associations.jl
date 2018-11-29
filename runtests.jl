if lowercase(get(ENV, "CI", "false")) == "true"
    include("install_dependencies.jl")
end

using TransferEntropy
using Test
using Distances


n_realizations = 50

@testset "Visitation frequency estimator" begin
	include("transferentropy/test_transferentropy_visitfreq.jl")
end

@testset "kNN (Kraskov) estimator" begin
	include("transferentropy/test_transferentropy_kraskov.jl")
end

@testset "Transfer operator grid estimator" begin
	include("transferentropy/test_transferentropy_transferoperator_grid.jl")
end
