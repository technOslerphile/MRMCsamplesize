# MRMCsamplesize
**An R Package to Estimate Sample Sizes for Multi-Reader Multi-Case (MRMC) Studies**

*MRMCsamplesize* is an R package with a function to estimate number of cases and readers required for a planned MRMC study. It does so based on the Obuchowski-Rockette (OR) model for statistical analysis of MRMC data. The package
is meant only for sample size estimation and it cannot be used for statistical analysis of MRMC data.

[*Available now in CRAN*](https://cran.r-project.org/web/packages/MRMCsamplesize/index.html)

[*Detailed documentation now available as a preprint in medRxiv*](https://www.medrxiv.org/content/10.1101/2023.09.25.23296069v1)

NOTE: The methodology does not incorporate the changes proposed by [*Hillis et. al.*](https://onlinelibrary.wiley.com/doi/abs/10.1002/sim.2532) by using a modified denominator degrees of freedom. The current methodology uses the denominator degrees of freedom proposed by Obuchowski and Rockette in their original publication. This is known to be very conservative. Option to incorporate modified degrees of freedom by Hillis is in development.

**How to install the package**
```
#You can install it directly from CRAN as usual

install.packages("MRMCsamplesize")

#You can also install using devtools the version in Github (this version is the same as in CRAN)
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

Dennis Robert, Saigopal Sathyamurthy S, Preetham Putha. MRMCsamplesize: An R Package for Estimating Sample Sizes for Multi-Reader Multi-Case Studies. medRxiv. Published online 2023. doi:10.1101/2023.09.25.23296069



