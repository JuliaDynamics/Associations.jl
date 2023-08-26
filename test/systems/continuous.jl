using StateSpaceSets: StateSpaceSet
using DynamicalSystemsBase: ContinuousDynamicalSystem
@test ContinuousDefinition <: SystemDefinition

@test trajectory(system(ChuaCircuitsBidir6()), 10, Δt = 0.05)[1] isa StateSpaceSet{6}
@test trajectory(system(ChuaScrollSine3()), 10, Δt = 0.05)[1] isa StateSpaceSet{3}
@test trajectory(system(HindmarshRose3()), 10, Δt = 0.05)[1] isa StateSpaceSet{3}
@test trajectory(system(LorenzBidir6()), 10, Δt = 0.05)[1] isa StateSpaceSet{6}
@test trajectory(system(LorenzForced9()), 10, Δt = 0.05)[1] isa StateSpaceSet{9}
@test trajectory(system(LorenzTransitive9()), 10, Δt = 0.05)[1] isa StateSpaceSet{9}
@test trajectory(system(MediatedLink9()), 10, Δt = 0.05)[1] isa StateSpaceSet{9}
@test trajectory(system(Repressilator6()), 10, Δt = 0.05)[1] isa StateSpaceSet{6}
@test trajectory(system(RosslerBidir6()), 10, Δt = 0.05)[1] isa StateSpaceSet{6}
@test trajectory(system(RosslerForced9()), 10, Δt = 0.05)[1] isa StateSpaceSet{9}
@test trajectory(system(RosslerLorenzUnidir6()), 10, Δt = 0.05)[1] isa StateSpaceSet{6}
@test trajectory(system(Thomas3()), 10, Δt = 0.05)[1] isa StateSpaceSet{3}
