using CausalityTools
using Test 

# Common errors
def = CMIRenyiJizba()

# Can't decompose CMIRenyiJizba/TERenyiJizba into Shannon entropies.
@test_throws ArgumentError EntropyDecomposition(def, LeonenkoProzantoSavani(Shannon()))
@test_throws ArgumentError EntropyDecomposition(def, PlugIn(Shannon()), OrdinalPatterns(m=2))
