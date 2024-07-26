using Test
using CausalityTools
defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname=defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "CausalityTools.jl" begin
    include("test_systems.jl")
    testfile("deprecations.jl")
    testfile("methods/methods.jl")
    testfile("independence/independence.jl")
    testfile("causal_graphs/causal_graphs.jl")
end
