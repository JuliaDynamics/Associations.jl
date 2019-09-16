using Reexport

@reexport module IntegrationDynamicalSystems
    import ..CausalityTests
    import CausalityToolsBase: causality
    import CausalityToolsBase: CausalityTest
    import CausalityToolsBase
    import DynamicalSystems: DynamicalSystem, DiscreteDynamicalSystem, trajectory, ContinuousDynamicalSystem

    ############################
    # Dynamical system resampling
    ############################
    include("DynamicalSystemResampling.jl")

    ############################
    # Dynamical system setup
    ############################
    include("DynamicalSystemSetup.jl")
    include("DiscreteSystemSetup.jl")

    ##########################
    # Continuous system setup
    ##########################
    include("ContinuousSystemSetup.jl")

    ############################################
    # causality tests for discrete systems
    ############################################

    """ 
        causality(system::DiscreteDynamicalSystem, setup::DiscreteSystemSetup, 
            test::CausalityTest)

    Apply the causality `test` to the given discrete `system` using the provided `setup` parameters.

    ## Example 

    ```julia
    # Use one of the built-in dynamical systems that ship with CausalityTools, 
    # `logistic2_unidir`, which is a system of two chaotic logistic maps `x` and `y` with 
    # coupling `x` to `y`, with `c_xy` controlling the coupling strength.
    sys = logistic2_unidir(c_xy = 1.0)

    # We treat the first variable as the source and the second variable as the target.
    # Iterate the system until we have time series with 200 observations. 
    setup = DiscreteSystemSetup(source = 1, target = 2, n_pts = 200)

    # Define some causality tests (with default values, you'd want to tweak the 
    # parameters to your needs)
    test_ccm = ConvergentCrossMappingTest(timeseries_lengths = [45, 50], n_reps = 20)
    test_cm = CrossMappingTest(n_reps = 10)
    test_vf = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = 1:5)
    test_tog = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = 1:5)
    test_jdd = JointDistanceDistributionTest()
    test_jddt = JointDistanceDistributionTTest()
    predtest = VisitationFrequencyTest(binning = RectangularBinning(5), ηs = -5:5)
    test_pa = PredictiveAsymmetryTest(predictive_test = predtest)

    # Compute the causality statistics using the different tests
    causality(sys, setup, test_ccm)
    causality(sys, setup, test_cm)
    causality(sys, setup, test_vf)
    causality(sys, setup, test_tog)
    causality(sys, setup, test_jdd)
    causality(sys, setup, test_jddt)
    ```
    """
    function causality(system::DiscreteDynamicalSystem, setup::DiscreteSystemSetup{NoResampling}, 
            test::CausalityTest)

        # Iterate the system and create orbit
        dt, Ttr, n_pts = setup.dt, setup.Ttr, setup.n_pts
        orbit = trajectory(system, n_pts, dt = dt, Ttr = Ttr)

        # Causality analysis 
        source = orbit[:, setup.source]
        target = orbit[:, setup.target]

        # Resampling scheme is `NoResampling`, so call causality without any constraints.
        causality(source, target, test)
    end

    """ 
        causality(system::ContinuousDynamicalSystem, setup::ContinuousSystemSetup, 
            test::CausalityTest)

    Apply the causality `test` to the given continuous `system` using the provided `setup` parameters.

    ## Example 

    Let's re-create the `lorenz_lorenz` system that ships with `CausalityTools` and analyse the 
    dynamical influence between some of its components.

    ```julia
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

    # Predictive asymmetry causality test 
    predtest = TransferOperatorGridTest(binning = RectangularBinning(5), ηs = -5:5)
    test_pa = PredictiveAsymmetryTest(predictive_test = predtest)

    causality(sys, setup, test_pa)
    ```
    """
    function causality(system::ContinuousDynamicalSystem, setup::ContinuousSystemSetup{NoResampling},
            test::CausalityTest)

        dt, Ttr, n_pts, sample_step = setup.dt, setup.Ttr, setup.n_pts, setup.sample_step
        
        # the system is recorded at times t0:dt:T
        T = n_pts*dt*sample_step
        
        orbit = trajectory(system, T, dt = dt, Ttr = Ttr*dt)[1:sample_step:end-1, :]
        #(percent_noise > 0) ? add_noise(o, percent_noise) : o

        # Causality analysis 
        source = orbit[:, setup.source]
        target = orbit[:, setup.target]
 
        # Resampling scheme is `NoResampling`, so call causality without any constraints.
        causality(source, target, test)
    end

    export causality, ContinuousSystemSetup, DiscreteSystemSetup
end

"""
    IntegrationDynamicalSystems

Module providing an interface for causality detection for dynamical systems as defined 
in DynamicalSystems.jl.
"""
IntegrationDynamicalSystems
