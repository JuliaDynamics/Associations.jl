using Test
using CausalityTools
defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname=defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "CausalityTools.jl" begin
    testfile("core.jl")
    testfile("methods/methods.jl")
    testfile("utils.jl")
    testfile("independence/independence.jl")
    testfile("causal_graphs/oce.jl")
    testfile("systems/systems.jl")
end

include("complexity_measures/test_effort_to_compress.jl")
include("complexity_measures/test_dynamical_complexity.jl")
include("complexity_measures/test_icc.jl")

#include("integrations/test_uncertaindata_integration.jl")
