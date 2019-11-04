import UncertainData:
    AbstractUncertainIndexValueDataset,
    ConstrainedIndexValueResampling,
    SequentialResampling,
    SequentialInterpolatedResampling,
    constrain,
    resample 

function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest,
    resampling::ConstrainedIndexValueResampling{N_VARIABLES, N_DATASETS}) where {N_VARIABLES, N_DATASETS}
    @warn "Performing causality test on uncertain index-value datasets without explicitly providing sequential sampling constraints. Results may be nonsensical, because the values are not sampled according to increasing indices."
    
    N_DATASETS != 2 ? error("Too many constraints ($N_DATASETS) provided. This `causality` method only takes `source` and `target` as data inputs. Thus, you should provide an `ConstrainedIndexValueResampling` containing two (`index_constraints::CONSTRAINT_TYPES, value_constraints::CONSTRAINT_TYPES)`-tuples. For example, `ConstrainedIndexValueResampling((TruncateStd(0.5), TruncateQuantiles(0.3, 0.7)), (TruncateStd(0.7), TruncateStd(1.1)))` would work.") : nothing
    
    # Take the (index_constraints, values_constraints)-tuple for the source and target
    # and represent as separate `IndexValueResampling` instances.
    cs_src = ConstrainedIndexValueResampling(resampling[1])
    cs_targ = ConstrainedIndexValueResampling(resampling[2])
  
    constrained_source = constrain(source, cs_src)
    constrained_target = constrain(target, cs_targ)
    
    [causality(resample(constrained_source.values), resample(constrained_target.values), test) for i = 1:resampling.n]
end


function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest,
    constrained_resampling::ConstrainedIndexValueResampling{N_VARIABLES, N_DATASETS},
    sequential_resampling::SequentialInterpolatedResampling{S, G}) where {N_VARIABLES, N_DATASETS, S, G}
    
    N_DATASETS != 2 ? error("Too many constraints ($N_DATASETS) provided. This `causality` method only takes `source` and `target` as data inputs. Thus, you should provide an `ConstrainedIndexValueResampling` containing two (`index_constraints::CONSTRAINT_TYPES, value_constraints::CONSTRAINT_TYPES)`-tuples. For example, `ConstrainedIndexValueResampling((TruncateStd(0.5), TruncateQuantiles(0.3, 0.7)), (TruncateStd(0.7), TruncateStd(1.1)))` would work.") : nothing
    
    # Take the (index_constraints, values_constraints)-tuple for the source and target
    # and represent as separate `IndexValueResampling` instances.
    cs_src = ConstrainedIndexValueResampling(constrained_resampling[1])
    cs_targ = ConstrainedIndexValueResampling(constrained_resampling[2])
    
    constrained_source = constrain(source, cs_src)
    constrained_target = constrain(target, cs_targ)
    
    [causality(resample(constrained_source, sequential_resampling)[2], 
               resample(constrained_target, sequential_resampling)[2], 
               test) for i = 1:constrained_resampling.n]
end

function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest,
    constrained_resampling::ConstrainedIndexValueResampling{N_VARIABLES, N_DATASETS},
    sequential_resampling::SequentialResampling{S}) where {N_VARIABLES, N_DATASETS, S}
    
    N_DATASETS != 2 ? error("Too many constraints ($N_DATASETS) provided. This `causality` method only takes `source` and `target` as data inputs. Thus, you should provide an `ConstrainedIndexValueResampling` containing two (`index_constraints::CONSTRAINT_TYPES, value_constraints::CONSTRAINT_TYPES)`-tuples. For example, `ConstrainedIndexValueResampling((TruncateStd(0.5), TruncateQuantiles(0.3, 0.7)), (TruncateStd(0.7), TruncateStd(1.1)))` would work.") : nothing
    
    # Take the (index_constraints, values_constraints)-tuple for the source and target
    # and represent as separate `IndexValueResampling` instances.
    cs_src = ConstrainedIndexValueResampling(constrained_resampling[1])
    cs_targ = ConstrainedIndexValueResampling(constrained_resampling[2])
    
    constrained_source = constrain(source, cs_src)
    constrained_target = constrain(target, cs_targ)
    
    [causality(resample(constrained_source, sequential_resampling)[2], 
               resample(constrained_target, sequential_resampling)[2], 
               test) for i = 1:constrained_resampling.n]
end