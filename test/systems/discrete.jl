using StateSpaceSets: Dataset
using DynamicalSystemsBase: DiscreteDynamicalSystem

@test Henon3() isa DiscreteSystem
@test ChaoticMaps3() isa DiscreteSystem

@test system(Henon3()) isa DiscreteDynamicalSystem
@test system(ChaoticMaps3()) isa DiscreteDynamicalSystem

n = 50
@test trajectory(system(Henon3()), n) isa Dataset
@test trajectory(system(ChaoticMaps3()), n) isa Dataset
