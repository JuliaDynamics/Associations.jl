using CausalityTools
using Test
using Random
rng = MersenneTwister(1234)

x, y = rand(rng, 100), rand(rng, 100)

# ----------------------------------------
# `MIShannon`
# ----------------------------------------
# On vector-valued inputs
def = MIShannon()
@test information(KSG1(def, k = 2), x, y) isa Real
@test information(KSG2(def, k = 2), x, y) isa Real
@test information(GaoOhViswanath(def, k = 2), x, y) isa Real
@test information(GaoKannanOhViswanath(def, k = 2), x, y) isa Real
@test information(GaussianMI(def), x, y) isa Real

# On `StateSpaceSet`s
data = [rand(rng, 50, 2) for i = 1:2]
x, y = StateSpaceSet.(data)
def = MIShannon()
@test information(KSG1(def, k = 2), x, y) isa Real
@test information(KSG2(def, k = 2), x, y) isa Real
@test information(GaoOhViswanath(def, k = 2), x, y) isa Real
@test information(GaoKannanOhViswanath(def, k = 2), x, y) isa Real
@test information(GaussianMI(def), x, y) isa Real
