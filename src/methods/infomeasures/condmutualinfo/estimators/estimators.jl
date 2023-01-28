include("FPVP.jl")
include("Rahimzamani.jl")
include("PoczosSchneiderCMI.jl")
include("MesnerShalisi.jl")

# Definition is actually never used, but we need to define it, so that calling `estimate`
# within independence tests work.
estimate(definition, measure::CMIShannon, est::ConditionalMutualInformationEstimator,
    x, y, z) = estimate(measure, est, x, y, z)

#include("TsallisCMIFuruichi.jl")
