#' Descriptive Statistics
#'
#' Calculates the mean, median, maximum, minimum, standard deviation, skewness, kurtosis & Jarque-Bera Normality test
#' @param x a Nx1 vector
#' @return descriptive statistics & p-values
#' @export
#' @author Alexandros Gabrielsen 
dstats <-function(x) {
  
T=length(x)
meanv=mean(x)
medianv=median(x)
maxv=max(x)
minv=min(x)
sigma=sd(as.numeric(x))
skew=sum(((x-meanv)/sigma)^3)/T
kurt=sum(((x-meanv)/sigma)^4)/T

tmp.test = jarquebera(x)
ans = as.matrix(t(c(meanv, medianv, maxv, minv, sigma, skew, kurt, tmp.test[1], tmp.test[2])))
colnames(ans) = c("Mean", "Median", "Maximum", "Minimum", "Std. Dev.", "Skewness", "Kurtosis", "Jarque-Bera", "Probability")
rownames(ans)="Statistics"

return=ans
}
