# Invariant measures

## From precomputed transfer operators.
When the transfer operator has been computed for a state space discretization, we can
derive an invariant measure (probability density) over the discrete states of the system. 
There are some different ways of doing this from observational data.

To obtain an invariant measure from a precomputed transfer operator, use the 
[`invariantmeasure`](invariantmeasure.md) function.

## From scratch, storing all relevant information.
If you want to keep track of all the pieces of information that goes into the estimation 
of an invariant measure, you may use one of the following estimators.
Currently, the following estimators are implemented.

- [`rectangularinvariantmeasure`](rectangularinvariantmeasure.md). Computes the invariant measure at a single resolution.
- [`InducedInvariantMeasure`](inducedinvariantmeasure.md). Computes the invariant measure for a partition at a target resolution induced by another resolution.
- [`AverageInvariantMeasure`](averageinvariantmeasure.md). Computes the invariant measure for a partition at a target resolution induced by multiple other partitions. 


