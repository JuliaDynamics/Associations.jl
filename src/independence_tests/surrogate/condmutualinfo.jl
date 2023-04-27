function SurrogateTest(measure::ConditionalMutualInformation, est::Nothing, args...;
        kwargs...)
    T = typeof(measure)
    txt = "Estimator not provided for measure $T. Cannot construct `SurrogateTest`\n" *
        "Do e.g. `SurrogateTest(CMIShannon(), FPVP())`"
    throw(ArgumentError(txt))
end
