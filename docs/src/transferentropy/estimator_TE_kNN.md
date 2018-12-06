# Transfer entropy (TE) estimators

## k Nearest neighbours (kNN) estimator
The kNN estimator computes transfer entropy as the sum of two mutual information terms, which are computed using the Kraskov estimator of mutual information ([Kraskov et al., 2004](https://doi.org/10.1103/PhysRevE.69.066138)). Implemented
for [Diego et al. (2018)](https://arxiv.org/abs/1811.01677).

## Documentation
```@docs
tekNN
```

## References
Kraskov, A., Stögbauer, H., & Grassberger, P. (2004). Estimating mutual information. Physical Review E - Statistical Physics, Plasmas, Fluids, and Related Interdisciplinary Topics. [https://doi.org/10.1103/PhysRevE.69.066138](https://doi.org/10.1103/PhysRevE.69.066138).

Diego, D., Agasøster Haaga, K., & Hannisdal, B. (2018, November 1). Transfer entropy computation using the Perron-Frobenius operator. Eprint ArXiv:1811.01677. Retrieved from [https://arxiv.org/abs/1811.01677](https://arxiv.org/abs/1811.01677).
