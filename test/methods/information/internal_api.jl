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
