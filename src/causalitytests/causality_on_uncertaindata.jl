
import UncertainData: 
    resample,
    AbstractUncertainIndexValueDataset,
    SequentialSamplingConstraint,
    RegularGrid

function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest)
    
    @warn """
    You're running a causality test by resampling `source` and `target`, which both have 
    uncertain indices, without ensuring proper ordering of the indices. If the 
    distributions/populations furnishing the indices have overlapping supports, you are 
    not guaranteed to have the correct ordering, and the results you get might be meaningless! 
    
    Use the following method instead:
    
    `causality(source::AbstractUncertainIndexValueDataset, 
        target::AbstractUncertainIndexValueDataset, 
        constraint::SequentialSamplingConstraint, 
        grid::RegularGrid)`. 
    
    For example, `causality(source, target, StrictlyIncreasing(), RegularGrid(start, stop, step))`
    will first pick a strictly increasing draw of the age model, interpolate to a regular grid 
    with spacing `step`, then run the test on the data interpolated to this grid, and finally
    perform the causality test.
    """
    idxs_s, vals_s = resample(source)
    idxs_t, vals_t = resample(target)
    causality(vals_s, vals_t, test)
end

function causality(source::AbstractUncertainIndexValueDataset, 
    target::AbstractUncertainIndexValueDataset, 
    test::CausalityTest,
    constraint::SCT,
    grid::RegularGrid) where SCT <: SequentialSamplingConstraint
    
    idxs_s, vals_s = resample(source, constraint, grid)
    idxs_t, vals_t = resample(target, constraint, grid)
    causality(vals_s, vals_t, test)
end