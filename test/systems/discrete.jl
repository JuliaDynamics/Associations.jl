using StateSpaceSets: Dataset
using DynamicalSystemsBase: DiscreteDynamicalSystem
@test DiscreteDefinition <: SystemDefinition

@test Anishchenko() isa DiscreteDefinition
@test AR1Unidir() isa DiscreteDefinition
@test AR1Bidir() isa DiscreteDefinition
@test Henon2() isa DiscreteDefinition
@test Henon3() isa DiscreteDefinition
@test Ikeda2() isa DiscreteDefinition
@test ChaoticMaps3() isa DiscreteDefinition
@test LinearMap2() isa DiscreteDefinition
@test Logistic2Unidir() isa DiscreteDefinition
@test Logistic2Bidir() isa DiscreteDefinition
@test Logistic3CommonDriver() isa DiscreteDefinition
@test Logistic4Chain() isa DiscreteDefinition
@test Nonlinear3() isa DiscreteDefinition

@test system(Anishchenko()) isa DiscreteDynamicalSystem
@test system(AR1Unidir()) isa DiscreteDynamicalSystem
@test system(AR1Bidir()) isa DiscreteDynamicalSystem
@test system(Henon2()) isa DiscreteDynamicalSystem
@test system(Henon3()) isa DiscreteDynamicalSystem
@test system(Ikeda2()) isa DiscreteDynamicalSystem
@test system(ChaoticMaps3()) isa DiscreteDynamicalSystem
@test system(LinearMap2()) isa DiscreteDynamicalSystem
@test system(Logistic2Unidir()) isa DiscreteDynamicalSystem
@test system(Logistic2Bidir()) isa DiscreteDynamicalSystem
@test system(Logistic3CommonDriver()) isa DiscreteDynamicalSystem
@test system(Logistic4Chain()) isa DiscreteDynamicalSystem
@test system(Nonlinear3()) isa DiscreteDynamicalSystem

n = 50
@test trajectory(system(Anishchenko()), n) isa Dataset
@test trajectory(system(AR1Unidir()), n) isa Dataset
@test trajectory(system(AR1Bidir()), n) isa Dataset
@test trajectory(system(Henon2()), n) isa Dataset
@test trajectory(system(Henon3()), n) isa Dataset
@test trajectory(system(Ikeda2()), n) isa Dataset
@test trajectory(system(ChaoticMaps3()), n) isa Dataset
@test trajectory(system(LinearMap2()), n) isa Dataset
@test trajectory(system(Logistic2Unidir()), n) isa Dataset
@test trajectory(system(Logistic2Bidir()), n) isa Dataset
@test trajectory(system(Logistic3CommonDriver()), n) isa Dataset
@test trajectory(system(Logistic4Chain()), n) isa Dataset
@test trajectory(system(Nonlinear3()), n) isa Dataset
