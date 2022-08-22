
abstract type SlidingWindow <: CausalityEstimator end
export ConstantWidthSlidingWindow
"""
    ConstantWidthSlidingWindow(estimator; n = 20, step = 1) <: SlidingWindow

Apply `estimator` to chunks of the time series, using windows of size `n` and `step` points 
between each window.
"""
struct ConstantWidthSlidingWindow{E, J} <: SlidingWindow
    estimator::E
    width::J
    step::J

    function ConstantWidthSlidingWindow(estimator::E; width::J = 20, step::J = 1) where {E, J}
        new{E, J}(estimator, width, step)
    end
end

ConstantWidthSlidingWindow(1, width = 20, step = 1)