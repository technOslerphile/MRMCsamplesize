# MRMCsamplesize
**An R Package to Estimate Sample Sizes for Multi-Reader Multi-Case (MRMC) Studies**

[//]: # [**Click here to go to the online R Package Documentation**](https://technoslerphile.github.io/autoCovariateSelection/index.html)

*MRMCsamplesize* is an R package with a function to estimate number of cases and readers required for a planned MRMC study. It does so based on the Obuchowski-Rockette (OR) model for statistical analysis of MRMC data. The package
is meant only for sample size estimation and it cannot be used for statistical analysis of MRMC data.

CRAN submission in progress. This will be updated when this package becomes available from CRAN.

**How to install the package**
```
#using devtools, you can install it from GitHub directly
library("devtools") 
install_github("technOslerphile/MRMCsamplesize")
```
**Usage example**
```
library("MRMCsamplesize")
result1 <- sampleSize_MRMC(endpoint = 'auc',J = 10,delta = 0.10,theta = 0.75,
                           rangeb = 0.1, rangew = 0.05, R = 1, r1 = 0.47,corr = FALSE)
result2 <- sampleSize_MRMC(endpoint = 'auc',J = 20,delta = 0.05,theta = 0.75,
                           rangeb = 0.2, rangew = 0.05, R = 1, r1 = 0.47,corr = TRUE, ICC = 0.5, s = 1.25)
result3 <- sampleSize_MRMC(endpoint = 'se',J = 15, delta = 0.05, theta = 0.75,
                           rangeb = 0.2, rangew = 0.025, R = 1, r1 = 0.5, corr = TRUE, ICC = 0.5, s = 1.25)
```
--------------------------------------------------------------------------------------------------------
I gladly welcome all suggestions for improvements and collaborations

[*Please report any bugs or issues here*](https://github.com/technOslerphile/MRMCsamplesize/issues)

**Please cite this package if you use it for your research**

  Robert D (2023). MRMCsamplesize: Sample Sizes for Multi-Reader Multi-Case (MRMC) Studies. R package version 1.0.0,
  <https://github.com/technOslerphile/MRMCsamplesize>.

