include("KLDivergence.jl")
include("RenyiDivergence.jl")
include("HellingerDistance.jl")
include("VariationDistance.jl")

abstract type MutualInformation <: BivariateInformationMeasure end
include("MIShannon.jl")
include("MITsallisMartin.jl")
include("MITsallisFuruichi.jl")
include("MIRenyiSarbu.jl")
include("MIRenyiJizba.jl")

abstract type JointEntropy <: BivariateInformationMeasure end
 # q-logarithm for Tsallis and Renyi joint entropies
logq(x, q) = (x^(1-q) - 1) / (1 - q)

include("JointEntropyShannon.jl")
include("JointEntropyRenyi.jl")
include("JointEntropyTsallis.jl")