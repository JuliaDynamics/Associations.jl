###############
# Discrete maps
###############
@testset "Discrete example systems" begin 

    # Check that initialisation happens correctly.
    @test isa(ar1_unidir(), DiscreteDynamicalSystem)
    @test isa(ar1_bidir(), DiscreteDynamicalSystem)
    @test isa(anishchenko1(), DiscreteDynamicalSystem)
    @test isa(henon2(), DiscreteDynamicalSystem)
    @test isa(ikeda(), DiscreteDynamicalSystem)
    @test isa(nonlinear3d(), DiscreteDynamicalSystem)
    @test isa(ulam(2), DiscreteDynamicalSystem)
    @test isa(logistic2_unidir(), DiscreteDynamicalSystem)
    @test isa(logistic2_bidir(), DiscreteDynamicalSystem)
    @test isa(logistic3(), DiscreteDynamicalSystem)
    @test isa(logistic4(), DiscreteDynamicalSystem)
    @test isa(var1(), DiscreteDynamicalSystem)
    @test isa(verdes(), DiscreteDynamicalSystem)

    # Discrete maps with lags directly return a Dataset
    @test isa(nontrivial_pegiun(n = 10), Dataset)
    @test isa(henon_triple(n = 10), Dataset)

    # Initialise all the systems and generate trajectories
    @test isa(trajectory(ar1_unidir(), 10), Dataset)
    @test isa(trajectory(ar1_bidir(), 10), Dataset)

    @test isa(trajectory(anishchenko1(), 10), Dataset)
    @test isa(trajectory(henon2(), 10), Dataset)
    @test isa(trajectory(ikeda(), 10), Dataset)
    @test isa(trajectory(nonlinear3d(), 10), Dataset)
    @test isa(trajectory(ulam(2), 10), Dataset)
    @test isa(trajectory(logistic2_unidir(), 10), Dataset)
    @test isa(trajectory(logistic2_bidir(), 10), Dataset)

    @test isa(trajectory(logistic3(), 10), Dataset)
    @test isa(trajectory(logistic4(), 10), Dataset)
    @test isa(trajectory(var1(), 10), Dataset)
    @test isa(trajectory(verdes(), 10), Dataset)

end