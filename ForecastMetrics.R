ForecastsMetrics <- function(Actual,Forecasts) {
#  Mean Square Error, MSE
#  Mean Absolute Deviation, MAD
#  Mean Logarithm of Absolute Errors, MLAE
#  Heteroskedasticity-adjusted Mean Square Error, HMSE 
#  Heteroskedasticity-adjusted Mean Absolute Error, HMAE 
#  Median Absolute Error, MAE
#  Median Absolute Percentage Error, MAPE
#  R2LOG
#  QLIKE
#  Success Ratio, SR
#MAPE
#MdAPE 
#RMSPE 
#RMdSPE 

#Symmetric Percentage Forecasts Metrics (Makridakis,1993)
#  sMAPE 
#  sMdAPE 
# Alexandros Gabrielsen   
  
#We suppose that Actual and Forecasts has been filtered
# to match same tenors beforehand
n<- length(Actual)
MFE = sum(Actual - Forecasts)/n 
MAE = sum(abs(Actual - Forecasts))/n 
MSE = sum((Actual - Forecasts)^2)/n
RMSE = sqrt(MSE)
UStat = RMSE/(sqrt(sum(Actual^2)/n)  + sqrt(sum(Forecasts^2)/n)) 
PCSP = length(which(Actual*Forecasts > 0))/n
#R2LOG = mean((log(Actual^2/Forecasts))^2)
#QLIKE = mean(log(Forecasts) + Actual^2/Forecasts)

#DA directional accuracy computation
P = length(which(Actual > 0))
hat_P = length(which(Forecasts > 0))
SRI =P*hat_P + (1- P)*(1- hat_P)
DA = (PCSP - SRI)/(sqrt(var(PCSP) - var(SRI)))

#Percentage Forecasts Metrics
MAPE = mean(abs(Actual - Forecasts)/abs(Actual))*100
MdAPE = median(abs(Actual - Forecasts)/abs(Actual))*100
RMSPE = sqrt(mean((abs(Actual - Forecasts)/abs(Actual))^2)*100)
RMdSPE = sqrt(median((abs(Actual - Forecasts)/abs(Actual))^2)*100)

#Symmetric Percentage Forecasts Metrics (Makridakis,1993)
sMAPE =  mean(abs(Actual-Forecasts)/(Actual+Forecasts)) *200
sMdAPE = median(abs(Actual-Forecasts)/(Actual+Forecasts)) *200

Labels <- c("MFE","MAE","MSE","RMSE","UStat","PCSP","MAPE","MdAPE","RMSPE","RMdSPE","sMAPE", "sMdAPE")
Values <- c(MFE,MAE,MSE,RMSE,UStat,PCSP,MAPE,MdAPE,RMSPE,RMdSPE,sMAPE,sMdAPE)

return(cbind(Labels,Values))
}

ForecastsMetricsBenchMark <- function(Actual,Forecasts,BenchMark) {
#We suppose that Actual and Forecasts, BenchMark has been filtered
# to match same tenors beforehand

r <- (Actual - Forecasts)/(Actual - BenchMark)
MRAE  = mean(abs(r))
MdRAE = median(abs(r))
GMRAE =gm_mean(abs(r))

Labels <- c("MRAE","MdRAE","GMRAE")
Values<- c(MRAE,MdRAE,GMRAE)

return(cbind(Labels,Values))
}

#Hyndman => relative in-sample
ForecastsMetricsInSample <- function(Actual,Forecasts,InSample) {
#We suppose that Actual and Forecasts has been filtered
# to match same tenors beforehand

MAEInSample = sum(abs(InSample[2:length(InSample)] - InSample[1:length(InSample) -1]))/(length(InSample) -1)

Q = (Actual - Forecasts)/MAEInSample

MASE = mean(abs(Q))
MSSE = mean((Q)^2)
RMSSE = sqrt(MSSE)

Labels <- c("MASE","MSSE","RMSSE")
Values<- c(MASE,MSSE,RMSSE)

return(cbind(Labels,Values))
}


