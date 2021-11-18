CalcPower\_tutorial
================

<style type="text/css">
  body{
  font-size: 16pt;
}
</style>

## Introduction

This function uses simulation to perform power analysis. It is designed
to explore the power of biological experiments and to suggest an optimal
number of experimental variables with reasonable power. The backbone of
the function is based on simr package, which fits a fixed effect or
mixed effect model based on the observed data and simulates response
variables. Users can test the power of different combinations of
experimental variables and parameters.

## Installation

``` r
library(devtools)
install_github('gladstone-institutes/CalcPower', build_vignettes=TRUE)
```

#### First few lines of the input data

``` r
head(CalcPower_data1)
#>   experiment plate  line classification feature1 feature2   feature3 feature4
#> 1       exp7     1 line1              0       62 3261.238   16.69571 118.4444
#> 2       exp7     1 line1              0       77 2033.726   54.10482 114.0946
#> 3       exp7     1 line1              0       56 1935.731 -100.00000 122.7329
#> 4       exp7     1 line1              0       54 1533.812 -100.00000 118.7066
#> 5       exp7     1 line1              0       19 1437.420 -100.00000 117.0279
#> 6       exp7     1 line1              0       30 1899.334 -100.00000 117.9523
#>   feature5
#> 1       10
#> 2       14
#> 3       13
#> 4       13
#> 5       19
#> 6       11
```

## Ex1. Varaince estimation and power calculation from pilot data

##### Levels will be increased to max(observed level)x5. For example, there are 9 experiments in the original data:

exp1, exp2, exp3, …, exp9

##### and we will test 9 x 5 experiments

#### code example:

``` r
calculate_power(data=CalcPower_data1,power_curve=1,
                 variance_estimate_from="data",condition_variable="classification",
                 experimental_variable=c("experiment","plate","line"), response_variable="feature1",
                 nsimn=10, target_parameters="experiment", levels=1)  
```

##### \* power\_curve=1 : to get a power curve that calculates power for different levels of the target parameter

##### \* variance\_estimate\_from =“data” : to estimate variance values from the input data

##### \* condition\_variable : cell status variable (ex control/case)

##### \* experimental\_variable : variables related to the experimental design

##### \* response\_varaible: phenotype

##### \* nsimn=10 : 10 iterations for power calculation ( try at least 100 or 1000 to get high accuracy)

##### \* target\_parameters=“experiment” : to explore “experiment”

##### \* levels = 1: to explore different levels of the target parameter

![Ex1 result](/Users/mingyoungshin/git/CalcPower/vignettes/ex1.jpeg)

## Ex2. Varaince estimation and power calculation from pilot data

##### Sample sizes will be increased to max(observed sample size)x5. For example, when there are max 71 samples per cell line in the original data:

``` r
table(CalcPower_data1$line)
#> 
#> line1 line2 line3 line4 line5 line6 line7 line8 
#>    42    40    37    51    71    37    42    20
```

##### and we will test max 71 x 5 samples per cell line

#### code example:

``` r
calculate_power(data=CalcPower_data1,power_curve=1,
                 variance_estimate_from="data",condition_variable="classification",
                 experimental_variable=c("experiment","plate","line"), response_variable="feature1",
                 nsimn=10, target_parameters="line", levels=0)  
```

###### \* levels = 0 : to explore different sample sizes of the target parameter

![Ex2 result](/Users/mingyoungshin/git/CalcPower/vignettes/ex2.jpeg)

## Ex3. User determined level count and output file name

##### Varaince estimation and power calculation (for a single level size) from pilot data with user determined level count and output file name. We will do Ex1 with max level size 15

#### code example:\`

``` r
calculate_power(data=CalcPower_data1,power_curve=0,
                 variance_estimate_from="data",condition_variable="classification",
                 experimental_variable=c("experiment","plate","line"), response_variable="feature1",
                 nsimn=10, target_parameters="experiment", levels=1, max_size=15,output="test.txt")  
```

##### \* power\_curve = 0 : to get the power for a single level of the target parameter

##### \* max\_size = 15 : to calculate power for the case when the target parameter has 15 levels (levels=1 therefore it will test levels not sample sizes)

##### \* output : to assign a name to the output file

##### Result:

<p class="comment">

Power for predictor ‘condition\_variable’, (95% confidence interval):  
80.00% (44.39, 97.48)  
  
Test: Likelihood ratio  
  
Based on 10 simulations, (0 warnings, 0 errors)  
alpha = 0.05, nrow = 544  
   
Time elapsed: 0 h 0 m 2 s  
  
nb: result might be an observed power calculation  

</p>

## Ex4. User determined effect size

#### Varaince estimation and power calculation from pilot data with user determined effect size of the condition variable

#### code example:

``` r
calculate_power(data=CalcPower_data1,power_curve=0,
                 variance_estimate_from="data",condition_variable="classification",
                 experimental_variable=c("experiment","plate","line"), response_variable="feature1",
                 nsimn=10, target_parameters="experiment", levels=1, effect_size = c(10))  
```

##### \* effect\_size = 10 : to assign 10 to the effect size of condition\_variable when the reesponse variable is “feature1”

##### Result:

<p class="comment">

Power for predictor ‘condition\_variable’, (95% confidence interval):  
50.00% (18.71, 81.29)  
  
Test: Likelihood ratio  
  
Based on 10 simulations, (0 warnings, 0 errors)  
alpha = 0.05, nrow = 1700  
  
Time elapsed: 0 h 0 m 10 s

</p>

## Ex5. Test two target parameters

#### Varaince estimation and power calculation from pilot data for two target parameters: Increase the number of levels of the first target parmeter and sample sizes of the second target parameter. We will test max 9 experiments and max 142 samples per cell line. There will be two results for each parameter

#### code example:

``` r
calculate_power(data=CalcPower_data1,power_curve=1,
                 variance_estimate_from="data",condition_variable="classification",
                 experimental_variable=c("experiment","plate","line"), response_variable="feature2",
                 nsimn=10, target_parameters=c("experiment","line"), levels=c(1,0), max_size=c(9,142))  
```

![Ex5-1
result](/Users/mingyoungshin/git/CalcPower/vignettes/ex5-1.jpeg)  
![Ex5-2 result](/Users/mingyoungshin/git/CalcPower/vignettes/ex5-2.jpeg)

## Ex6. Data with a single experimental category

#### If the pilot data only has a single experimental category for ‘experiment’,‘plat’ or’line’, we use known ICC values from other data.

#### code example:

Check sample size table of experimental varaibles:

``` r
table(CalcPower_data2$experiment,CalcPower_data2$plate,CalcPower_data2$line)
#> , ,  = line1
#> 
#>       
#>         1
#>   exp1 14
```

``` r
calculate_power(data=CalcPower_data2,power_curve=1,
                 variance_estimate_from="ICC",condition_variable="classification",
                 experimental_variable=c("experiment","plate","line"), response_variable="feature2",
                 nsimn=10, target_parameters=c("experiment"), levels=1, ICC=c(0.2,0.15,0.3))  
```

![Ex8 result](/Users/mingyoungshin/git/CalcPower/vignettes/ex8.jpeg)
