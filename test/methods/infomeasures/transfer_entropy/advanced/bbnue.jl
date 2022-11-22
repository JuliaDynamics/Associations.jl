import CausalityTools: construct_candidate_variables
using Entropies: VisitationFrequency, RectangularBinning

@testset "BBNUE" begin
     # Use periodic signals, so we also can test variable selection methods,
    # which for sensible testing, need to have their autocorrelation minima
    # > 1.
    s = sin.(1:100) .+ rand(100)
    t = sin.(1:100) .+ rand(100)
    c = sin.(1:100) .+ rand(100)

    # Internals that must work in order for the method to work.
    @testset "Variable exclusion" begin
        τexclude = 1
        vars = construct_candidate_variables([s], [t], τexclude = nothing)
        vars_ex = construct_candidate_variables([s], [t], τexclude = τexclude)
        fvars = Iterators.flatten(vars[1][1:end-1]) |> collect
        fvars_ex = Iterators.flatten(vars_ex[1][1:end-1]) |> collect
        @test τexclude ∈ abs.(fvars)
        @test τexclude ∉ abs.(fvars_ex)
        @test length(fvars) > length(fvars_ex)

        vars = construct_candidate_variables([s], [t], [c], τexclude = nothing)
        vars_ex = construct_candidate_variables([s], [t], [c], τexclude = τexclude)
        fvars = Iterators.flatten(vars[1][1:end-1]) |> collect
        fvars_ex = Iterators.flatten(vars_ex[1][1:end-1]) |> collect
        @test τexclude ∈ abs.(fvars)
        @test τexclude ∉ abs.(fvars_ex)
        @test length(fvars) > length(fvars_ex)
    end

    est = VisitationFrequency(RectangularBinning(3))
    te_st, params_st = bbnue(s, t, est, include_instantaneous = true)
    te_stc, params_stc = bbnue(s, t, c, est, include_instantaneous = false)
    @test te_stc isa Float64
    @test te_stc >= 0.0
end
