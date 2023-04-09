function SurrogateTest(measure::PMI, est::Nothing, args...;
        kwargs...)
    T = typeof(measure)
    txt = "Estimator not provided for measure $T. Cannot construct `SurrogateTest`\n" *
        "Do e.g. `SurrogateTest(PMI(), MesnerShalizi())`"
    throw(ArgumentError(txt))
end
