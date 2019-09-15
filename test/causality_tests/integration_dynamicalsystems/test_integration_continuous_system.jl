function eom_lorenz_lorenz(u, p, t)
    c_xy, c_yx, a₁, a₂, a₃, b₁, b₂, b₃ = (p...,)
    x1, x2, x3, y1, y2, y3 = (u...,)
    
    dx1 = -a₁*(x1 - x2) + c_yx*(y1 - x1)
    dx2 = -x1*x3 + a₂*x1 - x2
    dx3 = x1*x2 - a₃*x3
    dy1 = -b₁*(y1 - y2) + c_xy*(x1 - y1)
    dy2 = -y1*y3 + b₂*y1 - y2
    dy3 = y1*y2 - b₃*y3
    
    return SVector{6}(dx1, dx2, dx3, dy1, dy2, dy3)
end

function lorenz_lorenz(; u0 = rand(6), 
        c_xy = 0.2, c_yx = 0.0, 
        a₁ = 10, a₂ = 28, a₃ = 8/3, 
        b₁ = 10, b₂ = 28, b₃ = 9/3)
    ContinuousDynamicalSystem(eom_lorenz_lorenz, u0, [c_xy, c_yx, a₁, a₂, a₃, b₁, b₂, b₃])
end 

# Create an instance of the system with default parameters 
sys = lorenz_lorenz()

# Analysis setup. We'll check for an influence from variable x1 to y1, which are the 
# 1st and 4st components of the orbit. We'll generate orbits consisting of 500 points.
setup = ContinuousSystemSetup(source = 1, target = 4, n_pts = 500)

test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 20)
test_cm = CrossMappingTest(n_reps = 10)
test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)
test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1:5)
test_jdd = JointDistanceDistributionTest()
test_jddt = JointDistanceDistributionTTest()
predtest = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = -5:5)
test_pa = PredictiveAsymmetryTest(predictive_test = predtest)

# Run test 
@test causality(sys, setup, test_ccm) isa Vector{Vector{T}} where {T <: Real}
@test causality(sys, setup, test_cm) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_vf) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_tog) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_jdd) isa Vector{T} where {T <: Real}
@test causality(sys, setup, test_jddt) isa OneSampleTTest
