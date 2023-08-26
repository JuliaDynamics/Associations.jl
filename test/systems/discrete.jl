using StateSpaceSets: StateSpaceSet
using DynamicalSystemsBase: DiscreteDynamicalSystem
@test DiscreteDefinition <: SystemDefinition

@test Anishchenko() isa DiscreteDefinition
@test AR1Unidir() isa DiscreteDefinition
@test AR1Bidir() isa DiscreteDefinition
@test Henon2() isa DiscreteDefinition
@test Henon3() isa LaggedDiscreteDefinition
@test Ikeda2() isa DiscreteDefinition
@test ChaoticMaps3() isa DiscreteDefinition
@test ChaoticNoisyLinear2() isa LaggedDiscreteDefinition
@test Logistic2Unidir() isa DiscreteDefinition
@test Logistic2Bidir() isa DiscreteDefinition
@test Logistic3CommonDriver() isa DiscreteDefinition
@test Logistic4Chain() isa DiscreteDefinition
@test Nonlinear3() isa DiscreteDefinition

n = 50
@test trajectory(system(Anishchenko()), n)[1] isa StateSpaceSet
@test trajectory(system(AR1Unidir()), n)[1] isa StateSpaceSet
@test trajectory(system(AR1Bidir()), n)[1] isa StateSpaceSet
@test trajectory(system(Henon2()), n)[1] isa StateSpaceSet
@test trajectory(system(Henon3()), n)[1] isa StateSpaceSet
@test trajectory(system(Ikeda2()), n)[1] isa StateSpaceSet
@test trajectory(system(ChaoticMaps3()), n)[1] isa StateSpaceSet
@test trajectory(system(ChaoticNoisyLinear2()), n)[1] isa StateSpaceSet
@test trajectory(system(Logistic2Unidir()), n)[1] isa StateSpaceSet
@test trajectory(system(Logistic2Bidir()), n)[1] isa StateSpaceSet
@test trajectory(system(Logistic3CommonDriver()), n)[1] isa StateSpaceSet
@test trajectory(system(Logistic4Chain()), n)[1] isa StateSpaceSet
@test trajectory(system(Nonlinear3()), n)[1] isa StateSpaceSet
@test trajectory(system(Peguin2()), n)[1] isa StateSpaceSet
@test trajectory(system(UlamLattice()), n)[1] isa StateSpaceSet
@test trajectory(system(Var1()), n)[1] isa StateSpaceSet
@test trajectory(system(Verdes3()), n)[1] isa StateSpaceSet
