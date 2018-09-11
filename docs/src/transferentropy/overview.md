# Which transfer entropy estimators to choose?

For short time series, information flow estimators based on the transfer operator are the best bet. These compute transition probabilities in a more
rigorous way than visitation frequency based estimators. However, the transfer operator based approaches are computationally demanding. Exactly which estimator
to choose depends on the application:

- If you have a lot of data points (n > 1000), choose an estimator that computes transition probabilities based on counting visitation frequencies, or a nearest neighbor-based estimator. Estimators: `tetofreq` (recommended) and `tefreq`.
- For a moderate amount of points (100 < n < 1000), choose an estimator which computes transition probabilities from a triangulation-derived transfer operator using the *approximate simplex intersection approach*. The approximate simplex intersection is much faster than the exact, and generally doesn't deviate more than around 10 percent from the exact approach. Estimators: `tetotri`.
- For a moderate amount of points (100 < n < 1000), choose an estimator which computes transition probabilities from a triangulation-derived transfer operator using the *exact simplex intersection approach*. This method is the slowest, but yields the most accurate results you can get, given the information that is present in your data to begin with.
