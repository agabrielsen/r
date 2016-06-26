gof.reg <- function(resids, Y, p, LLF) {
  # PURPOSE:
  # Estimate the goodness of fit on residuals 
  #
  # INPUTS:
  # resids		a matrix of residuals
  # Y	        a matrix of the realised data
  # p         total number of explanatory variables, excluding intercept
  # LLF       LogLikelihood Function
  #
  # OUTPUTS:
  # a matrix covering a series of statistics such as
  # Residual Descriptive Statistics
  # Jarque-Bera
  # R^2
  # Adjusted R^2 
  # AIC      Akaike Information Criterion
  # BIC      Bayesian Information Criterion
  # SSE      Explained Sum of Squares
  # SST      Total Sum of Squares
  # SSR      Total Sum of Residuals
  # MSE      Mean Squared Error  
  # RMSE     Root Mean Square Error
  # ADF      Augmented-Dickey Fuller Unit Root Test
  # PP       Phillips-Perron Unit Root Test
  # KPSS     Kwiatkowski-Phillips-Schmidt-Shin Unit Root Test
  # LRTest   Likelihood Ratio Test
  # KS       Kolmogorov-Smirnov
  # SW       Shapiro-Wilk
  # CM       Crame Von Mises
  # AD       Anderson-Darling
  # SF       Shapiro-Francia
  # LB       Ljung-Box
  # BP       Box-Pierce
  #
  # Author
  # Alexandros Gabrielsen
  # 
  
  
  (R2 =  1 - sum(resids^2)/ sum((Y-mean(Y))^2))
  (R2Adj = 1-(1-R2)*(NROW(output[,X])-1)/(NROW(output[,X])-p))
  t.AIC = 2*c(p+2)-2*(as.numeric(LLF))
  t.BIC = -2*(LLF) +c(p+2)*log(NROW(Y))
  aa=t(dstats(resids))
  
  SSE =  sum(resids^2) # Explained Sum of Squares
  SST = sum((Y - mean(Y))^2)  # Total Sum of Squares
  SSR = sum((Y-resids -mean(Y))^2)
  
  #This statistic is also known as the fit standard error and the standard error of the regression. 
  #It is an estimate of the standard deviation of the random component in the data, and is defined as
  MSE = SSE/c(NROW(Y)-3) 
  RMSE = sqrt(MSE)
  
  # Augmented Dickey-Fuller
  ADF = adf.test(resids,k=1)
  
  # Philips-Perron
  PP = pp.test(resids)
  
  #KPSS
  KPSS = kpss.test(resids)
  
  # Kolmogorov-Smirnov Test
  KS = ks.test(resids,"pnorm",mean=0,sd=1) 
  
  # Likelihood Ratio Test
  # Fitting an intercept model
  bb = tryCatch(arima(Y, c(0,0,0),method="ML", optim.method="BFGS",optim.control = list(maxit = 2500, reltol=1e-8)), error=function( err ) FALSE, warning=function( warn ) FALSE ) # This checks if an error is provided
  if( !is.logical( bb ) ) {
    return = bb
  } else { return =  NULL }
  
  LRTest = 2*(LLF-bb$loglik)
  pval=NULL
  for (i in seq(0.001,0.999,0.001)) {
    pval = rbind(pval, c(1-i,qchisq(i, df=c(3))))
  }
  
  if (length(which(pval[,2]>LRTest))>0) {
    LRpval = pval[min(which(pval[,2]>LRTest)),1]
  } else {
    LRpval = 0.00001
  }
  rm(pval)
  
  # Cramer Von Mises
  CM = cvm.test(resids)
  
  # Anderson-Darling Test
  AD = ad.test(resids)
  
  # Shapiro-Wilk Test
  SW = shapiro.test(resids)
  
  # Shapiro-Francia
  SF = sf.test(resids)
  
  xx=NULL
  xx=rbind(xx,c("Mean",aa[1],"","",""))
  xx=rbind(xx,c("Median",aa[2],"","",""))
  xx=rbind(xx,c("Maximum",aa[3],"","",""))
  xx=rbind(xx,c("Minimum",aa[4],"","",""))
  xx=rbind(xx,c("Std. Dev.",aa[5],"","",""))
  xx=rbind(xx,c("Skewness",aa[6],"","",""))
  xx=rbind(xx,c("Kurtosis",aa[7],"","",""))
  xx=rbind(xx,c("Jarque-Bera",aa[8],"","",aa[9]))
  xx=rbind(xx,c("SSE",SSE,"","",""))            
  xx=rbind(xx,c("SST",SST,"","","")) 
  xx=rbind(xx,c("SSR",SSR,"","",""))
  xx=rbind(xx,c("MSE",MSE,"","",""))
  xx=rbind(xx,c("RMSE",RMSE,"","","")) 
  xx=rbind(xx,c("AIC",t.AIC,"","",""))
  xx=rbind(xx,c("BIC",t.BIC,"","",""))
  xx=rbind(xx,c("LogLike",LLF,"","",""))
  xx=rbind(xx,c("R-Squared",R2,"","",""))
  xx=rbind(xx,c("R-Squared Adj",R2Adj,"","",""))
  xx=rbind(xx,c("Augmented-Dickey Fuller",as.numeric(ADF$statistic),"","",as.numeric(ADF$p.value)))
  xx=rbind(xx,c("Phillips-Perron",as.numeric(PP$statistic),"","",as.numeric(PP$p.value)))
  xx=rbind(xx,c("Kwiatkowski-Phillips-Schmidt-Shin",as.numeric(KPSS$statistic),"","",as.numeric(KPSS$p.value)))
  xx=rbind(xx,c("Kolmogorov-Smirnov",as.numeric(KS$statistic),"","",as.numeric(KS$p.value)))
  xx=rbind(xx,c("Likelihood Ratio",as.numeric(LRTest),"","",as.numeric(LRpval)))
  xx=rbind(xx,c("Shapiro-Wilk",as.numeric(SW$statistic),"","",as.numeric(SW$p.value)))
  xx=rbind(xx,c("Crame Von Mises", as.numeric(CM$statistic),"","",as.numeric(CM$p.value)))
  xx=rbind(xx,c("Anderson-Darling",  as.numeric(AD$statistic),"","",as.numeric(AD$p.value)))
  xx=rbind(xx,c("Shapiro-Francia",  as.numeric(SF$statistic),"","",as.numeric(SF$p.value)))
  
  # Ljung-Box
  for (k in c(2):20) {
    bb = Box.test(resids, lag = k, type = c("Ljung-Box"), fitdf = c(k-1))
    xx=rbind(xx,c(paste("Ljung-Box",sep=""),bb$statistic,k,c(k-1), bb$p.value))
  }
  
  # Box-Pierce
  for (k in c(2):20) {
    bb = Box.test(resids, lag = k, type = c("Box-Pierce"), fitdf = c(k-1))
    xx=rbind(xx,c(paste("Box-Pierce",sep=""),bb$statistic,k,c(k-1), bb$p.value))
  }
  
  colnames(xx) = c("StatName", "StatValue", "StatLag", "StatDF", "StatPValue")
  
  cc = c("StatValue", "StatLag", "StatDF", "StatPValue")
 
  return(as.data.frame(xx))
}
