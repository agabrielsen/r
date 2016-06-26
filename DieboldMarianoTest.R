DieboldMarianoTest <- function(Actual,Forecast1, Forecast2, power_)
{
	loss1 <- Actual - Forecast1
	loss2 <- Actual - Forecast2
	testres <- dm.test(loss1, loss2, power= power_)
	
	Labels <- c("Statistic", "pvalue")
	Values <- c(testres$statistic, testres$p.value)
	
	return(cbind(Labels,Values))
}

