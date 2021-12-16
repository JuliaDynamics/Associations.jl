using Reexport 

@reexport module PredictiveAsymmetry 
    include("predictive_asymmetry.jl")
    include("automated_pa.jl")
    include("symmetric.jl")
end