using Test
using CausalityTools
defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname=defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "CausalityTools.jl" begin
    testfile("methods/methods.jl")
    testfile("utils.jl")
    testfile("independence/independence.jl")
    testfile("systems/systems.jl")
end

#include("integrations/test_uncertaindata_integration.jl")
