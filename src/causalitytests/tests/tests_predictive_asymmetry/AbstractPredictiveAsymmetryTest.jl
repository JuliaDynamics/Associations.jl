"""
    AbstractPredictiveAsymmetryTest

An abstract type for predictive asymmetry causality tests.

Concrete subtypes must implement the following fields:

- **`predictive_test`**. A predictive test to use for computing the 
    asymmetries. 
"""
abstract type AbstractPredictiveAsymmetryTest <: CausalityTest end


""" 
    get_ηs
    
Get the prediction lags for a predictive asymmetry test.
"""
get_ηs(x::CT) where {CT <: AbstractPredictiveAsymmetryTest} = x.predictive_test.ηs

export AbstractPredictiveAsymmetryTest, get_ηs

