
################################################################
# Distance based causality tests
################################################################
include("tests_distance_based/DistanceBasedCausalityTest.jl")

# S-measure test
# ---------------------
import ..s_measure

include("tests_distance_based/SMeasureTest.jl")

# Joint distances tests
# ---------------------
import ..joint_distance_distribution

include("tests_distance_based/JointDistancesCausalityTest.jl")
include("tests_distance_based/JointDistanceDistributionTest.jl")
include("tests_distance_based/JointDistanceDistributionTTest.jl")

# Cross mapping tests
# ---------------------
include("tests_distance_based/CrossMappingTest.jl")
include("tests_distance_based/ConvergentCrossMappingTest.jl")


################################################################
# Entropy based causality tests
################################################################
include("tests_entropy_based/EntropyBasedCausalityTest.jl")

# Transfer entropy causality tests
# ---------------------------------------
include("tests_entropy_based/TransferEntropyCausalityTest.jl")
include("tests_entropy_based/VisitationFrequencyTest.jl")
include("tests_entropy_based/NearestNeighbourMITest.jl")
include("tests_entropy_based/TransferOperatorGridTest.jl")
include("tests_entropy_based/ApproximateSimplexIntersectionTest.jl")
include("tests_entropy_based/ExactSimplexIntersectionTest.jl")

################################################################
# Predictive asymmetry causality tests
################################################################
include("tests_predictive_asymmetry/AbstractPredictiveAsymmetryTest.jl")
include("tests_predictive_asymmetry/PredictiveAsymmetryTest.jl")
include("tests_predictive_asymmetry/NormalisedPredictiveAsymmetryTest.jl")


# For the metatypes, we need to know the return types of the tests. These 
# are defined in the following file (should probably encode this as 
# a type parameter instead).
include("get_return_types.jl")
