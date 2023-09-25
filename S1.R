#######################################################################################################
# R code to:
#          1. Generate the examples in the main manuscript 
#          2. R output sample size numbers from Tables 1, 2 and 3 of supplementary material (S1)
#######################################################################################################

if(length(setdiff("MRMCsamplesize", rownames(installed.packages()))) > 0 ){
  install.packages("MRMCsamplesize")
  library("MRMCsamplesize")
}else{
  library("MRMCsamplesize")
}

##############################################################################
###                    EXAMPLE 1                                       ######
##############################################################################
ex1 <- sampleSize_MRMC(endpoint = 'auc',
                       J = 20, 
                       delta = 0.05 ,
                       theta = 0.75, 
                       r1 =0.47, 
                       rangeb = 0.20, 
                       rangew = 0.05)
ex1$ORSampleSizeResults
##############################################################################





##############################################################################
###                    EXAMPLE 2                                        ######
##############################################################################
ex2 <- sampleSize_MRMC(endpoint = 'auc',
                       J = 20, 
                       delta = 0.05 ,
                       theta = 0.75, 
                       r1 =0.47, 
                       rangeb = 0.20, 
                       rangew = 0.05,
                       corr = TRUE,
                       ICC = 0.5,
                       s = 1.25)
print(ex2$ORSampleSizeResults)
##############################################################################


##############################################################################
###                    Supplementary Tables 1, 2 and 3                  ######
##############################################################################

#Reference Sample Size Numbers in the Supplementary Tables are sourced from https://www.ajronline.org/doi/full/10.2214/ajr.175.3.1750603
#The R code only generates the R output numbers. The comparison of R output numbers with reference sample size numbers was done manually.

small <- c(0.01, 0.005)
moderate <- c(0.05, 0.025)
large <- c(0.1, 0.05)
auc <- c(0.75, 0.9)
delta <- c(0.05, 0.10, 0.15)
R <- c(1,2,4)
J <- c(4,6,10)

tab <- data.frame("TableName" = as.character(),
                  "J" = as.numeric(),
                  "AUC" = as.numeric(),
                  "delta" = as.numeric(),
                  "Ratio" = as.numeric(),
                  "variability" = as.character(),
                  "R_Output_Sample_Size" = as.numeric())

ntotal.small <- as.numeric()
ntotal.moderate <- as.numeric()
ntotal.large <- as.numeric()

iter <- 1
for(l in 1:length(J)){
  for (i in 1:length(auc)){
    for (j in 1:length(delta)){
      for(k in 1:length(R)){
        
        res.small <- tryCatch(sampleSize_MRMC(endpoint = 'auc',
                                              J = J[l],
                                              delta = delta[j],
                                              theta = auc[i],
                                              r1 =0.47,
                                              rangeb = small[1],
                                              rangew = small[2],
                                              R = R[k]),
                              error = function(e) NA)
        
        res.moderate <- tryCatch(sampleSize_MRMC(endpoint = 'auc',
                                                 J = J[l],
                                                 delta = delta[j],
                                                 theta = auc[i],
                                                 r1 =0.47,
                                                 rangeb = moderate[1],
                                                 rangew = moderate[2],
                                                 R = R[k]),
                                 error = function(e) NA)
        res.large <- tryCatch(sampleSize_MRMC(endpoint = 'auc',
                                              J = J[l],
                                              delta = delta[j],
                                              theta = auc[i],
                                              r1 =0.47,
                                              rangeb = large[1],
                                              rangew = large[2],
                                              R = R[k]),
                              error = function(e) NA)
        ntotal.small[iter] <- ifelse(!is.na(res.small), res.small$ORSampleSizeResults$nTotal, NA)
        ntotal.moderate[iter] <- ifelse(!is.na(res.moderate),res.moderate$ORSampleSizeResults$nTotal, NA)
        ntotal.large[iter] <- ifelse(!is.na(res.large),res.large$ORSampleSizeResults$nTotal, NA)
        
        iter <- iter + 1
      }
    }
  }
}

R <- rep(rep(R, 6),3)
delta <- rep(c(rep(delta[1],3), rep(delta[2],3), rep(delta[3],3), rep(delta[1],3), rep(delta[2],3), rep(delta[3],3)),3)
auc <- rep(c(rep(auc[1],9), rep(auc[2],9)),3)
J <- c(rep(J[1],18), rep(J[2], 18), rep(J[3],18))
tab.name <- c(rep("Table 1",18), rep("Table 2", 18), rep("Table 3",18))


table.validation <- data.frame("Table Name" = tab.name,
                               "J" = J,
                               "AUC" = auc,
                               "delta" = delta,
                               "Ratio" = R,
                               "R_Output_Sample_Size_Small_Reader_Variability)`" = ntotal.small,
                               "R_Output_Sample_Size_Moderate_Reader_Variability)`" = ntotal.moderate,
                               "R_Output_Sample_Size_Large_Reader_Variability)`" = ntotal.large)

