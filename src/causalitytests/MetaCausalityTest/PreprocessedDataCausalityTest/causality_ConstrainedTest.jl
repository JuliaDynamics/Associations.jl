"""
    causality(source::AbstractUncertainIndexValueDataset, 
        target::AbstractUncertainIndexValueDataset, 
        test::ConstrainedTest{CT, CR}) where {CT, CR}

Apply a causality test of type `CT` to `test.n` independent 
realisations of `source` and `target`, after first constraining 
the supports of the uncertain values furnishing the datasets.

See also [`ConstrainedTest`](@ref).
"""
function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::ConstrainedTest{CT, CR}) where {CT, CR}
    
     # Take the (index_constraints, values_constraints)-tuple for the source and target
    # and represent as separate `IndexValueResampling` instances.
    cs_src = ConstrainedIndexValueResampling(test.constrained_resampling[1])
    cs_targ = ConstrainedIndexValueResampling(test.constrained_resampling[2])

    constrained_source = constrain(source, cs_src)
    constrained_target = constrain(target, cs_targ)
    
    [causality(
            resample(constrained_source)[2], 
            resample(constrained_target)[2], 
            test.test) for i = 1:test.n]
end

export causality