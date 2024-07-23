using Test 
using CausalityTools

# ----------------------------------------------------------------
# This file tests internal functions.
# ----------------------------------------------------------------
def_renyi = CMIRenyiSarbu(; q = 5, base = 5)
def_tsallis = CMITsallisPapapetrou(; q = 5, base = 5)
def_shannon = CMIShannon(; base = 5)
est_renyi = PlugIn(Renyi(; q = 0.5, base = 2))
est_tsallis = PlugIn(Tsallis(; q = 0.5, base = 2))
est_shannon = PlugIn(Shannon(; base = 2))

new_est_renyi = CausalityTools.estimator_with_overridden_parameters(def_renyi, est_renyi)
new_est_tsallis = CausalityTools.estimator_with_overridden_parameters(def_tsallis, est_tsallis)
new_est_shannon = CausalityTools.estimator_with_overridden_parameters(def_shannon, est_shannon)
@test new_est_renyi == PlugIn(Renyi(; q = 5, base = 5)) 
@test new_est_tsallis == PlugIn(Tsallis(; q = 5, base = 5)) 
@test new_est_shannon == PlugIn(Shannon(; base = 5)) 


p1 = Probabilities([0.1, 0.2, 0.3])
p2 = Probabilities([0.1, 0.2, 0.3, 0.4])
@test_throws DimensionMismatch CausalityTools.size_match(KLDivergence(), p1, p2)