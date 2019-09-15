
##################################################################################################
# Integration with UncertainData.jl (sampling from full supports of the furnishing distributions)
#################################################################################################

import HypothesisTests: OneSampleTTest

# Vectors of uncertain values
uvals_x = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100]
uvals_y = [UncertainValue(Normal, rand(Normal(0, 5)), abs(rand(Normal(0, 3)))) for i = 1:100];

# UncertainValueDataset
UVX = UncertainValueDataset(uvals_x)
UVY = UncertainValueDataset(uvals_y)

# UncertainIndexDataset
UVX_idx = UncertainIndexDataset(uvals_x)
UVY_idx = UncertainIndexDataset(uvals_y)

# Real-valued vectors
x = resample.(uvals_x);
y = resample.(uvals_y);

jtest = JointDistanceDistributionTTest()

@test causality(x, y, jtest) isa OneSampleTTest
@test causality(uvals_x, uvals_y, jtest) isa OneSampleTTest
@test causality(x, uvals_y, jtest) isa OneSampleTTest
@test causality(uvals_x, y, jtest) isa OneSampleTTest
@test causality(UVX, UVY, jtest) isa OneSampleTTest
@test causality(x, UVY, jtest) isa OneSampleTTest
@test causality(UVX, y, jtest) isa OneSampleTTest


################################################################
# Integration with UncertainData.jl (with sampling constraints)
################################################################
onevar_constraints = ConstrainedResampling(TruncateStd(1))
twovar_constraints = ConstrainedResampling(TruncateStd(2), TruncateStd(1))

test_jddt = JointDistanceDistributionTTest() 
@test causality(x, y, twovar_constraints, test_jddt) isa OneSampleTTest 
@test causality(uvals_x, uvals_y, twovar_constraints, test_jddt) isa OneSampleTTest 
@test causality(x, uvals_y, onevar_constraints, test_jddt) isa OneSampleTTest 
@test causality(uvals_x, y, onevar_constraints, test_jddt) isa OneSampleTTest 
@test causality(UVX, UVY, twovar_constraints, test_jddt) isa OneSampleTTest 
@test causality(uvals_x, UVY, twovar_constraints, test_jddt) isa OneSampleTTest 
@test causality(UVX, uvals_y, twovar_constraints, test_jddt) isa OneSampleTTest 
@test causality(x, UVY, onevar_constraints, test_jddt) isa OneSampleTTest 
@test causality(UVX, y, onevar_constraints, test_jddt) isa OneSampleTTest 