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