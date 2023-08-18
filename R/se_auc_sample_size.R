#' Estimate sample sizes for standalone studies with sensitivity or AUC as endpoint
#' @description
#' \code{sampleSize_MRMC} This function returns number of cases required for a standalone study for endpoints of sensitivity and AUC.
#'
#' @author Dennis Robert \email{dennis.robert.nm@gmail.com}
#'
#' @param endpoint Character string to inform what is the endpoint (Figure-Of-Merit - FOM) of the standalone study. Values can be either \code{auc} or \code{sensitivity}.
#' @param theta Expected average value of the FOM Must be a value between 0 and 1.
#' @param precision Required precision of the point estimate of FOM. This is equivalent to half-width of the confidence interval. Must be a numeric value between 0 and 1.
#' @param R Ratio of non-diseased cases to diseased cases. Defaults to 1.
#' @param power Power to detect \code{delta} given all other assumptions. Default value is 0.8 corresponding to 80 percent power.
#' @param alpha The type I error rate. Default value is 0.05 corresponding to 5 percent type I error (significance level).
#' @param var_auc Variance estimation method when endpoint is \code{auc}. Defaults to the string \code{obuchowski}. If value is changed to \code{blume}, then method proposed by Blume (2009) will be used to estimate the variance.
#' @param corr Logical value indicating if \code{ICC (intra-cluster correlation)} has to be adjusted (\code{TRUE}) or not (\code{FALSE}). Defaults to \code{FALSE}.
#' @param ICC A numerical value between 0 and 1 indicating the expected ICC if \code{corr} is \code{TRUE}.
#' @param s Average number of lesions in diseased cases. This must be a numeric value greater than or equal to 1.
#' @return
#' A named list
#' \itemize{
#'   \item \code{SampleSizeResults} - A list containing the sample size results.
#'}
#'
#' @references
#' \itemize{
#' \item Flahault A, Cadilhac M, Thomas G. Sample size calculation should be performed for design accuracy in diagnostic test studies. J Clin Epidemiol. 2005 Aug;58(8):859-62. doi: 10.1016/j.jclinepi.2004.12.009. PMID: 16018921.
#' \item Zhou, X.-H., Obuchowski, N.A. and McClish, D.K. (2011). Sample Size Calculations. In Statistical Methods in Diagnostic Medicine (eds X.-H. Zhou, N.A. Obuchowski and D.K. McClish). https://doi.org/10.1002/9780470906514.ch6
#' }
#' @details
#' When \code{corr = FALSE}, the \code{nUnits_i} in \code{SampleSizeResults} is the number of diseased cases. The number of total cases (\code{nTotal}) required will depend on the
#' the ratio \code{R} specified.
#' When \code{corr = TRUE}, the anticipated correlation between units within the same diseased cases are adjusted and the \code{nUnits_i} in \code{SampleSizeResults}
#' list is the number of units in diseased cases assuming independence. The number of diseased cases required in this scenario will be given
#' by \code{nCases_c}. Again, \code{nTotal} will depend on the \code{R} specified.
#' @examples
#' library("MRMCsamplesize")
#' result1 <- sampleSize_Standalone(endpoint = "auc", theta = 0.9, precision = 0.05,
#'                                  R = 1, corr = TRUE, ICC = 0.5, s = 1.25)
#' result2 <- sampleSize_Standalone(endpoint = "Se", theta = 0.8, precision = 0.05, R = 1)
#' @export

sampleSize_Standalone <- function(endpoint = "auc",
                                  theta,
                                  precision,
                                  R = 1,
                                  power = 0.8,
                                  alpha = 0.05,
                                  var_auc = "obuchowski",
                                  corr = FALSE,
                                  ICC = NULL,
                                  s = NULL){


  if (precision > 0.20 & precision < 1) {warning("Precision seems to be very minimal. Are you sure to power the study for such a low precision of your estimate?")}
  if (theta >= 1 | theta <=0) {stop("The conjectured 'theta' value does not seem probable. It has to be a numeric value 0 < theta < 1")}
  if (precision > 1 | precision <= 0) {stop("Improbable precision. It has to be numeric value 0 < precision < 1")}



  if(tolower(endpoint) %in% c("se", "sens", "sensitivity") | var_auc == "blume"){
    var.theta= theta*(1-theta) #variance calculation for Se
  }else if(tolower(endpoint) %in% c("auc") & var_auc == "obuchowski"){
    A <-  stats::qnorm(theta)*1.414
    var.theta <- ((0.0099 * exp(-A^2/2)) * ( (5*A^2 +8) + (A^2 +8)/R))
  }else{
    stop("Endpoint has to be either 'sensitivity' or 'auc'")
  }

  num = (stats::qnorm(1-alpha/2) + stats::qnorm(power)) * sqrt(var.theta)
  den = precision
  n = ceiling((num/den)^2)
  nControls <- ceiling(n*R)

  n.total <- n + nControls

  if (!corr){
    text <- "No adjustment for ICC"
    nUnits_i <- n
    nCases_c <- NA
    nControls <- nControls
    nTotal <- n.total
    de <- NA
    s <- NA
  } else {
    if(is.null(ICC) | is.null(s)) {stop("ICC and s must be valid numeric values to take correlation within diseased cases into account")}
    text <- "Adjusted for ICC"
    s <- s
    de = 1 + (s-1)*ICC
    nUnits_i <- n
    nUnits_c = ceiling(nUnits_i * de)
    nCases_c = ceiling(nUnits_c/s)
    nControls = ceiling(nCases_c*R)
    nTotal = nCases_c + nControls
  }



  METHOD <- "Sample Size Estimation Results"

  NOTE <- paste0("\n",
                  "ICC: Is intra-class correlation (ICC) considered while estimating sample size", "\n",
                  "nUnits_i: Number of required units with presence of at least one lesion assuming independence between units", "\n",
                  "nCases_c: Number of required cases with presence of at least one lesion after adjusting for ICC", "\n",
                  "nControls: Number of required controls (images/patients without any lesion)", "\n",
                  "nTotal: Total sample size (cases)", "\n",
                  "DE: Design effect due to ICC", "\n",
                  "s: Assumed average number of lesions in cases with lesions", "\n",
                  "power: Assumed power", "\n",
                  "alpha: Significance level (Type I error rate)", "\n"
  )

  out <- structure(list(ICC = text,
                         nUnits_i = nUnits_i,
                         nCases_c = nCases_c,
                         nControls = nControls,
                         nTotal = nTotal,
                         DE = de,
                         s = s,
                         power = power,
                         alpha = alpha,
                         method = METHOD, note = NOTE), class = "power.htest")

  return(list("SampleSizeResults" = out))

}
