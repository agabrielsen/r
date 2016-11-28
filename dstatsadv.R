#' Descriptive Statistics, Normality & Univariate Tests
#'
#' Calculates the mean, median, maximum, minimum, standard deviation, skewness, kurtosis & Jarque-Bera Normality test, Kolmogorov-Smirnov,
#' Cramer Von Mises, Anderon-Darling, Shapiro-Wilk, Shapiro-Francia, ADF Unit Root, Phillips-Perron, KPSS
#' @param x a Nx1 vector
#' @return descriptive statistics & p-values
#' @export
#' @author Alexandros Gabrielsen 
dstatsadv <-function(x) {
  
T=length(x)
meanv=mean(x)
medianv=median(x)
maxv=max(x)
minv=min(x)
sigma=sd(as.numeric(x))
skew=sum(((x-meanv)/sigma)^3)/T
kurt=sum(((x-meanv)/sigma)^4)/T
tmp.test = jarquebera(x)

KS = ks.test(x,"pnorm",mean=mean(x),sd=sd(x))    # Kolmogorov-Smirnov Test
CM = cvm.test(x)                                 # Cramer Von Mises
AD = ad.test(x)                                  # Anderson-Darling Test
SW = shapiro.test(x)                             # Shapiro-Wilk Test
#lillie.test(zz)                                 # Lilliefors (Kolmogorov-Smirnov) test
SF = sf.test(x)                                  # Shapiro-Francia
adf= adf.test(x,k=1)                             # ADF Unit Root Testing
pp = pp.test(x)                                  # Phillips-Perron Unit Root
kpss = kpss.test(x)                              # KPSS Unit Root

ans = as.matrix(t(c(meanv, medianv, maxv, minv, sigma, skew, kurt, tmp.test[1], tmp.test[2],
                    as.numeric(KS$statistic),as.numeric(KS$p.value),
                    as.numeric(CM$statistic),as.numeric(CM$p.value),
                    as.numeric(AD$statistic),as.numeric(AD$p.value),
                    as.numeric(SW$statistic),as.numeric(KS$p.value), 
                    as.numeric(SF$statistic),as.numeric(SF$p.value),
                    adf$statistic,adf$p.value,
                    pp$statistic, pp$p.value,
                    kpss$statistic,kpss$p.value
                    )))
colnames(ans) = c("Mean", "Median", "Maximum", "Minimum", "Std. Dev.", "Skewness", "Kurtosis", 
                  "Jarque-Bera", "JBProb", "Kolmogorov-Smirnov", "KSProb",
                  "Cramer Von Mises", "CMProb", "Anderson-Darling", "ADProb",
                  "Shapiro-Wilk", "SWProb", "Shapiro-Francia", "SFProb", "ADF", "ADFProb",
                  "PP", "PPProb", "KPSS", "KPSSProb")
rownames(ans)="Statistics"

return=ans
}
