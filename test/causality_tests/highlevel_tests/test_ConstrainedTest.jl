
cx = (TruncateStd(1), TruncateStd(1))
cy = (TruncateStd(1), TruncateStd(1))
cs = ConstrainedIndexValueResampling(cx, cy)

# If n is not provided, set it to 1
ctest = ConstrainedTest(CrossMappingTest(), cs)
@test ctest isa ConstrainedTest
@test ctest.n == 1

ctest = ConstrainedTest(CrossMappingTest(), cs, 10)
@test ctest isa ConstrainedTest
@test ctest.n == 10