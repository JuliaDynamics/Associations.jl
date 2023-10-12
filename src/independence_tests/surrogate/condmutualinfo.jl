function SurrogateAssociationTest(measure::ConditionalMutualInformation, est::Nothing, args...;
        kwargs...)
    T = typeof(measure)
    txt = "Estimator not provided for measure $T. Cannot construct `SurrogateAssociationTest`\n" *
        "Do e.g. `SurrogateAssociationTest(CMIShannon(), FPVP())`"
    throw(ArgumentError(txt))
end
