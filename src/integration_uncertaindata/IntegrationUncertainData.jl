using Reexport

@reexport module IntegrationUncertainData
    import ..CausalityTests
    import CausalityToolsBase: causality
    import UncertainData: resample, AbstractUncertainValue, AbstractUncertainValueDataset, ConstrainedResampling
    import CausalityToolsBase: CausalityTest#, causality 
    
    include("integration_with_uncertaindata.jl")
    include("integration_with_uncertaindata_constrained_resampling.jl")
    
    # We may specialize behaviour if needed, but defining `causality`
    # for the generic tests works well for now, so we don't need to 
    # import these files.
    #include("integration_uncertaindata_ConvergentCrossMappingTest.jl")
    #include("integration_uncertaindata_CrossMappingTest.jl")
    #include("integration_uncertaindata_VisitationFrequencyTest.jl")
    #include("integration_uncertaindata_TransferOperatorGridTest.jl")
    #include("integration_uncertaindata_JointDistanceDistributionTest.jl")
    #include("integration_uncertaindata_JointDistanceDistributionTTest.jl")

    export causality

end