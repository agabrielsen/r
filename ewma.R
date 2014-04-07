# Exponentially Weighted Moving Average
# ----------------------------------------------------
# PURPOSE:
# Estimates the Exponentially Weighted Moving Average.
#---------------------------------------------------
# USAGE:
# [M, V, H, C] = ewma(data,lambda, method)
#
# INPUTS:
# data = ( M x N ) vector
# lambda = a scalar
# method = 0 to estimate univariate time series, 1 to estimate
# multivariate
#
# OUTPUTS:
# M = ( M x N ) mean vector
# V = ( M x N ) volatility vector
# H = ( M x M x N ) covariance vector
# C = ( M x M x N ) correlation vector
# R = ( M x M x N ) squared residuals
#
# M(t) = lamda*M(t-1) + (1-lamda)*X(t)
# H(t) = lamda*((data(t-1)-M(t-1)).^2) + (1-lamda)*H(t-1)
#
# NOTES:
# 1. the data vector is allowed to have different time steps
#---------------------------------------------------
# Author:
# Alexandros Gabrielsen
# Date: 04/2011
# v.1.0
#---------------------------------------------------

# EMWA Function
ewma = function(data, lambda, method) {

# Size of data vector
T = NROW(data)
N = NCOL(data)

# Prespecify vectors and variables
startdata = array(0, N) # Will store where data starts
M = array(NA, c(T, N)) # Mean vector
V = array(NA, c(T, N)) # Volatility vector
if (method == 1){
H = array(NA, c(N, N, T)) # Covariance vector
C = array(NA, c(N, N, T)) # Correlation vector
R = array(NA, c(N, N, T)) # Squared residuals
H1 = as.matrix(read.csv(paste(filepath,"Ht.csv", sep=""), header=FALSE))
H[ , , 1] = H1
C[ , , 1]= cov2cor(H1)
}

# Find all NA and NAN in the vector
W = is.na(data)

# Find where data starts
for (k in (1:N)) {
ww=data.frame(table(W[,k]))
startdata[k]=ww[2,2]
if (is.na(startdata[k])) {
startdata[k] = 0
}
}
startdata = startdata+1; # adjusting to start from correct time step

# Correct initialization and estimation
for (i in (1:N)){
# Initial values for the mean process
M[(startdata[i]),i] = mean(data[,i], na.rm = TRUE)
V[(startdata[i]),i] = var(data[,i], na.rm = TRUE)

# Estimating the mean process
for (j in ((startdata[i]+1):T)){
M[j,i]=(1-lambda)*data[j,i]+lambda*M[j-1,i]
V[j,i]=(1-lambda)*((data[j,i]-mean(data[,i], na.rm = TRUE))^2)+lambda*V[j-1,i]
}

if (method == 1){
# Correcting the covariance vector at time step 1
if (startdata[i] != 1) {
H[, i, 1] = c(NA)
}
}
}
rm(i, j, k ,ww, W) # clear unused variables

mm = colMeans(data, na.rm = TRUE)
# Estimating covariance and correlation vectors
if (method == 1){
for (i in (2:T)){
R[ , , i] = (data[i,]-mm)%*%t(data[i,]-mm)
H[ , , i] = (1-lambda)*R[ , ,i] + lambda*H[ , ,i-1]

# Correcting covariance vector at each time step
for (j in (1:N)){
if (startdata[j] == i) {
H[,j,i] = (1-lambda)*R[,j,i] + lambda*H1[,j]
H[j,,i] = t(H[,j,i])
}
}
# Estimating correlation vector
C[ , , i] = H[, , i] / (sqrt(diag(H[, , i]))%*%t(sqrt(diag(H[, , i]))))
}
rm(i,j) # clear unused variables
}

if (method == 1){
return=list(M, V, H, C)} else {
return=list(M, V)}

} # function ends
