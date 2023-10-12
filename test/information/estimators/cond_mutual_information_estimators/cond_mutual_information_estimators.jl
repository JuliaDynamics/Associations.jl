using CausalityTools
using Test
using Random
rng = MersenneTwister(1234)

x, y, z = rand(rng, 100), rand(rng, 100), rand(rng, 100);

# ----------------------------------------
# `CMIShannon`
# ----------------------------------------
# On vector-valued inputs
def = CMIShannon()
@test information(MesnerShalizi(def, k = 2), x, y, z) isa Real
@test information(FPVP(def, k = 2), x, y, z) isa Real
@test information(Rahimzamani(def, k = 2), x, y, z) isa Real
@test information(GaussianCMI(def), x, y, z) isa Real

# On `StateSpaceSet`s
data = [rand(rng, 50, 2) for i = 1:3]
x, y, z = StateSpaceSet.(data)
def = CMIShannon()
@test information(MesnerShalizi(def, k = 2), x, y, z) isa Real
@test information(FPVP(def, k = 2), x, y, z) isa Real
@test information(Rahimzamani(def, k = 2), x, y, z) isa Real
@test information(GaussianCMI(def), x, y, z) isa Real

# ----------------------------------------
# `CMIRenyiPoczos`
# ----------------------------------------

def = CMIRenyiPoczos()
@test information(PoczosSchneiderCMI(def, k = 2), x, y, z) isa Real

data = [rand(rng, 50, 2) for i = 1:3]
x, y, z = StateSpaceSet.(data)
@test information(PoczosSchneiderCMI(def, k = 2), x, y, z) isa Real
