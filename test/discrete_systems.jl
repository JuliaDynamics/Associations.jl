###############
# Discrete maps
###############

# Check that initialisation happens correctly.
@test isa(ar1(), DiscreteDynamicalSystem)
@test isa(anishchenko1(), DiscreteDynamicalSystem)
@test isa(henon2(), DiscreteDynamicalSystem)
@test isa(henon4(), DiscreteDynamicalSystem)
@test isa(linear3d_nonlinearcoupling(), DiscreteDynamicalSystem)
@test isa(logistic2(), DiscreteDynamicalSystem)
@test isa(logistic3(), DiscreteDynamicalSystem)
@test isa(logistic4(), DiscreteDynamicalSystem)
@test isa(var1(), DiscreteDynamicalSystem)
@test isa(var1coupled(), DiscreteDynamicalSystem)
@test isa(verdes(), DiscreteDynamicalSystem)

# Discrete maps with lags directly return a Dataset
@test isa(nontrivial_pegiun(n = 10), Dataset)
@test isa(henon_triple(n = 10), Dataset)

# Initialise all the systems and generate trajectories
@test isa(trajectory(ar1(), 10), Dataset)
@test isa(trajectory(anishchenko1(), 10), Dataset)
@test isa(trajectory(henon2(), 10), Dataset)
@test isa(trajectory(henon4(), 10), Dataset)
@test isa(trajectory(linear3d_nonlinearcoupling(), 10), Dataset)
@test isa(trajectory(logistic2(), 10), Dataset)
@test isa(trajectory(logistic3(), 10), Dataset)
@test isa(trajectory(logistic4(), 10), Dataset)
@test isa(trajectory(var1(), 10), Dataset)
@test isa(trajectory(var1coupled(), 10), Dataset)
@test isa(trajectory(verdes(), 10), Dataset)
