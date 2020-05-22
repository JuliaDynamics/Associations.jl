# Invariantizing a set of points

In this package, unless otherwise specified, invariance refers to the situation when a set
of points have been discretized by a triangulation.

## Forward linear map of triangulation simplices
A triangulation consists a set of disjoint simplices filling up the volume inside the
convex hull of the points. Some algorithms in this package[^1] first triangulates a set
of points, then use the forward linear map of the simplices (move each simplex vertex one step forward in time according to the observed orbit) to do some calculation.

These algorithms triangulate the first ``N-1`` points ``p_1, p_2, \ldots, p_{N-1}``, leaving
the last point untouched. This way, every vertex ``p_i`` are mapped to ``p_{i + 1}``,
including ``p_{N-1}`` which is mapped to ``p_{N}`` (this would not have worked if we triangulated all points).


### Ensuring that there is no information loss
For these algorithms to yield correct results, all simplices must be mapped to simplices
that lie within the convex hull of ``p_1, p_2, \ldots, p_{N-1}``. If Should ``p_{N}`` lies
outside the convex hull of the preceding points, information is lost when estimating
transition probabilities. Thus, we need a way of making sure that is the case. The
[`forwardlinearmap_invariant`](../../glossary/invariantizing/forwardlinearmap_invariant.md) function can be used to check
this condition.


## Invariantizing a set of points
There are several ways of making sure ``p_{N}`` lies inside the convex hull of the preceding
points. If your dataset is large enough, it might be sufficient to just remove points
from the end of the dataset until invariance is achieved.

Another way is to manipulate the last point by moving it towards the center of the convex
hull of the points. This will, of course, introduce some bias, but for very sparse and
noisy data, this might be a better option than removing points from an already limited
dataset. The [`invariantize`](invariantize.md) function does precisely this.


[^1]:
    Specifically, the [exact](../../transferoperator/transferoperator_triang_exact.md) and [approximate](../../transferoperator/transferoperator_triang_approx.md) triangulation
    transfer operator estimators, and the corresponding transfer entropy estimators.
