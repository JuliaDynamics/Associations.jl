# Examples of coupled dynamical systems

## Testing the performance of causality detection methods
The algorithms in this library are intended to measure the information flow, or dynamical coupling strength, between time series. To test the performance of the algorithms, we can apply them on computed generated time series for which we know the equations of motion.
We proceed as follows.

1. First, pretend we don't know the structure and parameters of the underlying equations from which we have generated time series.
2. Then, by looking only at the time series and the dynamical information intrinsic to the time series, can we retrieve the directionality of interactions the underlying system?

Numerous examples of coupled dynamical systems appear in the causality detection literature. Here, we provide examples of discrete and continuous coupled systems used to test causality algorithms.

### Other available system to study
*Note that the `DynamicalSystems` package also implements examples of dynamical systems that have been studied intensively in the literature. However, not many of these systems are coupled in a manner that is well-suited for causality detection algorithms. This requires dynamical systems that are coupled such that we can control the relative dynamical influence between variables.*
