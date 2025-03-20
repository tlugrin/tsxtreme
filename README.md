# tsxtreme R package

## A New Approach to Fitting Time Series Extremal Dependence

Classical approaches to fitting extremes of time series with short-term
dependence use a pre-processing stage where independent extremes are filtered
out of the original series. A well-known approach is the peaks-over-threshold
method which typically involves the definition of arbitrary clusters of
exceedances over a high threshold. The distribution fitted to the cluster maxima
depends on the exact definition of the underlying clusters and can yield
badly-biased estimates of related risk measures.

This package implements a Bayesian semiparametric approach (Lugrin et al., 2016,
_Extremes_) based on the conditional tail approach of Heffernan and Tawn (2004,
_JRSSB_) and formalised by Heffernan and Resnick (2007, _Ann. Appl. Probab._).
For comparison purposes, the package offers additional methods to compute
risk measures using an empirical and a stepwise approach using maximum
likelihood.

## Documentation

More details and a gentle introduction to the functionality of the package can
be found in the `vignette("time-series-extremes")` of the package.

## Installation

```
# Install release version from CRAN}
install.packages("tsxtreme")

# Install development version from GitHub
devtools::install("tlugrin/tsxtreme")
```

## References

A non-exhaustive list of references to papers whose results are used and
implemented in this package.

* Lugrin, T., Davison, A. C. and Tawn, J. A. (2016) Bayesian uncertainty
  management in temporal dependence of extremes. _Extremes_, **19**, 491--515.

* Coles, S., Heffernan, J. E. and Tawn, J. A. (1999) Dependence measures for
  extreme value analyses. _Extremes_, **2**, 339--365.
  
* Keef, C., Papastathopoulos, I. and Tawn, J. A. (2013) Estimation of the
  conditional distribution of a multivariate variable given that one of its
  components is large: Additional constraints for the Heffernan and Tawn model.
  _Journal of Multivariate Analysis_, **115**, 396--404.

* Heffernan, J. E. and Tawn, J. A. (2004) A conditional approach for
  multivariate extreme values. _Journal of the Royal Statistical Society
  Series B_, **66**, 497--546.

* Heffernan, J. E. and Resnick, S. I. (2007) Limit laws for random vectors
  with an extreme component. _Annals of Applied Probability_, **17**,
  537--571.

* Ledford, W. A. and Tawn, J. A. (2003) Diagnostics for dependence within time
  series extremes. _Journal of the Royal Statistical Society Series B_,
  **65**, 521--543.

* Davison, A. C. and Smith, R. L. (1990) Models for exceedances over high
  thresholds. _Journal of the Royal Statistical Society Series B_, **52**,
  393--442.
