# Available transfer entropy estimators


## Naive visitation frequency based estimator
The simplest transfer entropy estimator (`transferentropy_visitfreq`, or its
alias `tefreq`) assumes that the long-term probabilities of occupying states
can be estimated by the visitation frequency of each bin in the partitioned
state space. Visitation frequencies are determined by counting.

```@docs
transferentropy_visitfreq
```

## Transfer operator based estimators

The visitation frequency based estimator
(`transferentropy_transferoperator_visitfreq`, or its alias `tetofreq`) is the
fastest of the transfer operator transfer entropy estimators, and should be
used for any datasets of 1000 points or more. It works in a similar manner
as the `tefreq` estimator, but also takes into account the information about
the dynamics provided by the transfer operator.


```@docs
transferentropy_transferoperator_visitfreq
```

The triangulation based estimator (`transferentropy_transferoperator_triang`, or
its alias `tetotri`) is slower, but is computed based on a transfer operator
computed from a triangulation, which better captures the dynamics of the
underlying system.

```@docs
transferentropy_transferoperator_triang
```
