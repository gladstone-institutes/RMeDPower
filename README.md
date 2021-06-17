# calcPower

## How to install
```
install.packages("devtools")
library(devtools)
install_github('gladstone-institutes/calcPower', build_vignettes=TRUE)
library(calcPower)
```

## Functions

* `calculate_power` - This function uses simulation to perform power analysis. It is designed to explore the power of biological experiments and to suggest an optimal number of experimental variables with reasonable power. The backbone of the function is based on [simr](https://cran.r-project.org/web/packages/simr/index.html) package, which fits a fixed effect or mixed effect model based on the observed data and simulates response variables. Users can test the power of different combinations of experimental variables and parameters.
