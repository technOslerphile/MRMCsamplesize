% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrmc_sample_size.R
\name{sampleSize_MRMC}
\alias{sampleSize_MRMC}
\title{Estimate sample sizes for MRMC studies}
\usage{
sampleSize_MRMC(
  endpoint = "auc",
  J,
  delta,
  rangeb,
  rangew,
  theta,
  R = 1,
  r1,
  r2 = 0.3,
  r3 = 0.3,
  rb = 0.8,
  K = 1,
  power = 0.8,
  alpha = 0.05,
  nu1 = 1,
  var_auc = "obuchowski",
  reader_var_estimation_method = "normal",
  n_reading_sessions_per_reader = 2,
  corr = FALSE,
  ICC = NULL,
  s = NULL
)
}
\arguments{
\item{endpoint}{Character string to inform what is the Figure-of-Merit (FOM) which will be used as the endpoint of the MRMC study. Values can be either \code{auc} or \code{sensitivity}.}

\item{J}{The number of readers for the study. It is recommended to have minimum 5 readers in any MRMC study.}

\item{delta}{Effect size denoting the anticipated difference in the endpoint between the two interventions/imaging-modalities/techniques. Typically chosen values are 0.04, 0.05 and 0.06. Should be between 0 and 1.}

\item{rangeb}{Inter-reader variability range (sometimes referred to as between-reader variability) denoting the anticipated difference between the highest accuracy of any reader in the study and the lowest accuracy of any reader in the study. Should be a numeric value between 0 and 1.}

\item{rangew}{Intra-reader variability range (sometimes referred to as within-reader variability) denoting the anticipated difference between the accuracies of a reader who interprets the same images using the same imaging technique at two different times. Should be a numeric value between 0 and 1.}

\item{theta}{Expected average value of the FOM for the \code{J} readers.}

\item{R}{Ratio of non-diseased cases to diseased cases. Defaults to 1.}

\item{r1}{Correlation between FOMs of readers when same cases are evaluated by the same reader using different modalities.}

\item{r2}{Correlation between FOMs when the same cases are evaluated by different readers using the same modality. It is assumed that \code{r2 = r3} for default calculations.}

\item{r3}{Correlation between FOMs when the same cases are evaluated by different readers using different modalities. It is assumed that \code{r2 = r3} for default calculations.}

\item{rb}{Correlation between FOMs when the same readers evaluate cases using different modalities. The default value is \code{0.8}.}

\item{K}{Number of times each reader interprets the same case from the same modality. This is always equal to 1 in a fully-crossed paired-reader paired-case study design with two modalities.}

\item{power}{Power to detect \code{delta} given all other assumptions. Default value is 0.8 corresponding to 80 percent power.}

\item{alpha}{The type I error rate. Default value is 0.05 corresponding to 5 percent type I error (significance level).}

\item{nu1}{Numerator degrees of freedom of the F-distribution which will be used to estimate the non-centrality parameter (lambda).}

\item{var_auc}{Variance estimation method when endpoint is \code{auc}. Defaults to the string \code{obuchowski}. If value is changed to \code{blume}, then method proposed by Blume (2009) will be used to estimate the variance.}

\item{reader_var_estimation_method}{A value = \code{normal} uses the assumption that the accuracy of readers are distributed normally and thus the relationship between range and standard deviation can be used to estimate the inter and intra reader variances from \code{rangeb} and \code{rangew}. Any other value will use a rule of thumb to estimate inter and intra reader variances by dividing \code{rangeb} and \code{rangew} by 4 followed by squaring it. \code{normal} method is typically more conservative especially when J is less than 30-35.}

\item{n_reading_sessions_per_reader}{Number of times each reader interprets each case. Defaults to 2 which corresponds to a typical MRMC study with 2 modalities.}

\item{corr}{Logical value indicating if \code{ICC (intra-cluster correlation)} has to be adjusted (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{FALSE}.}

\item{ICC}{A numerical value between 0 and 1 indicating the expected ICC if \code{corr} is \code{TRUE}.}

\item{s}{Average number of lesions in diseased cases. This must be a numeric value greater than or equal to 1.}
}
\value{
A list within a list object with two named lists
\itemize{
  \item \code{varComponents} - A list containing the estimated values of the OR variances and correlation components.
  \item \code{ORSampleSizeResults} - A list containing the sample size results.
}
}
\description{
\code{sampleSize_MRMC} This function returns number of cases required for a MRMC study for a given number of readers.
}
\details{
When \code{corr = FALSE}, the \code{nUnits_i} in \code{ORSampleSizeResults} list is the number of diseased cases. The number of total cases (diseased + non-diseased; \code{nTotal}) required will depend on the
the ratio \code{R} specified.
When \code{corr = TRUE}, the anticipated correlation between units within the same diseased cases are adjusted and the \code{nUnits_i} in \code{ORSampleSizeResults}
list is the number of units in diseased cases assuming independence. The number of diseased cases required in this scenario will be given
by \code{nCases_c}. Again, \code{nTotal} required will depend on the \code{R} specified.
}
\examples{
library("MRMCsamplesize")
result1 <- sampleSize_MRMC(endpoint = 'auc',J = 10,delta = 0.10,theta = 0.75,
rangeb = 0.1, rangew = 0.05, R = 1, r1 = 0.47,corr = FALSE)
result2 <- sampleSize_MRMC(endpoint = 'auc',J = 20,delta = 0.05,theta = 0.75,
rangeb = 0.2, rangew = 0.05, R = 1, r1 = 0.47,corr = TRUE, ICC = 0.5, s = 1.25)
result3 <- sampleSize_MRMC(endpoint = 'se',J = 15, delta = 0.05, theta = 0.75,
rangeb = 0.2, rangew = 0.025, R = 1, r1 = 0.5, corr = TRUE, ICC = 0.5, s = 1.25)
}
\references{
\itemize{
\item Obuchowski NA, Hillis SL. Sample size tables for computer-aided detection studies. AJR Am J Roentgenol. 2011 Nov;197(5):W821-8. doi: 10.2214/AJR.11.6764. PMID: 22021528; PMCID: PMC3494304
\item Obuchowski NA. & Rockette HE. (1995) Hypothesis testing of diagnostic accuracy for multiple readers and multiple tests an anova approach with dependent observations, Communications in Statistics - Simulation and Computation, 24:2, 285-308, DOI: 10.1080/03610919508813243
\item Obuchowski NA. Sample size tables for receiver operating characteristic studies. AJR Am J Roentgenol. 2000;175(3):603-608. doi:10.2214/ajr.175.3.1750603
\item Rockette HE, Campbell WL, Britton CA, Holbert JM, King JL, Gur D. Empiric assessment of parameters that affect the design of multireader receiver operating characteristic studies. Acad Radiol. 1999;6(12):723-729. doi:10.1016/s1076-6332(99)80468-1
\item Blume JD. Bounding Sample Size Projections for the Area Under a ROC Curve. J Stat Plan Inference. 2009 Mar 1;139(1):711-721. doi: 10.1016/j.jspi.2007.09.015. PMID: 20160839; PMCID: PMC2631183.
}
}
\author{
Dennis Robert \email{dennis.robert.nm@gmail.com}
}
