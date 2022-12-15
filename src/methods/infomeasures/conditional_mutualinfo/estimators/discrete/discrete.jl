
struct ShannonCMI{E <: Renyi} <: ConditionalMutualInformation
    e::E
    function CMIShannon(e::E) where E <: Renyi
        e.q ≈ 1.0 || error("CMIShannon not defined Renyi entropy with q=$(e.q)")
        new{E}(e)
    end
end

function estimate(def::ShannonCMI, est::ProbabilitiesEstimator)
