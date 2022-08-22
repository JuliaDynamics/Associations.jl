""" The supertype of all time series causality estimators. """
abstract type CausalityEstimator end

include("utils/sliding_window.jl")