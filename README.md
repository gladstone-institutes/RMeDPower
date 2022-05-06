# RMeDPower

## How to install
```
install.packages("devtools")
library(devtools)
install_github('gladstone-institutes/RMeDPower', build_vignettes=TRUE)
library(RMeDPower)
```

## Functions

* `calculate_power` - This function uses simulation to perform power analysis. It is designed to explore the power of biological experiments and to suggest an optimal number of experimental variables with reasonable power. The backbone of the function is based on [simr](https://cran.r-project.org/web/packages/simr/index.html) package, which fits a fixed effect or mixed effect model based on the observed data and simulates response variables. Users can test the power of different combinations of experimental variables and parameters.

* 'calculate_lmer_estimates' - This function performs a linear mixed model analysis using lmer.

* 'transform_data' - This functions makes quantile-quanitle (qq) plots of i) raw residual values ii) log-transformed residual values iii) raw residual values after removing outliers, and iv) log-transformed residual values after removing outliers. To detect outliers, the function uses Rosner's test.

* 'check_normality' - This function makes a quantile-quantile (qq) plot of the residual values of the mixed effects model. Users can check the normality of residual values by examining qqplot.
