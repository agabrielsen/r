VaRLR <-function(fdata, VaR, alpha, position) { 
# 
#----------------------------------------------------------------------- 
# PURPOSE:  
# Value-at-Risk backtesting for long and short positions using 
# the unconditional coverage, independence and conditional coverage family  
# of tests 
#----------------------------------------------------------------------- 
# USAGE:  
# results = VaRLR(fdata, VaR, alpha, position, options) 
# 
# INPUTS: 
# fdata:     ( m x 1 ) vector of the out-of-sample data, 
# VaR:       ( m x 1 ) vector of VaR estimates 
# alpha:     a% 
# position:  Long or Short positions  
#----------------------------------------------------------------------- 
# OUTPUTS:  
# results:  PF:      Percentage of Failures 
#           TUFF:    Time Until First Failure 
#           LRTUFF:  Likelihood Ratio of Time Until First Failure 
#           LRUC:    Likelihood Ratio Unconditional Coverage   
#           LRIND:   Likelihood Ratio Independence Coverage 
#           LRCC:    Likelihood Ratio Conditional Coverage    
#           Basel:   Basel II Accord 
#----------------------------------------------------------------------- 
# REFERENCES: 
# BASEL II. (2005). "International convergence of capital measurement and 
#        captial standards.", Basel Committee on Banking Supervision. 
# Christoffersen, P., (1998). "Evaluating interval forecasts." 
#         International Economic Review, 39, 841-862. 
# Christoffersen, P., (2003). "Elements of Financial Risk Management."  
#          Academic Press, Elsevier Science 
# Kupiec. P., (1995). "Techniques for Verifying the Accuracy of Risk  
#          Management Models." Journal of Derivatives, 3,73.84. 
#----------------------------------------------------------------------- 
# Author: 
# Alexandros Gabrielsen 
# Date: 01/2014 
#----------------------------------------------------------------------- 

if (position == "Long") { hit = fdata<VaR } else { hit = fdata>VaR } 
n1 = sum(hit)         # Number of Violations 
n0 = NROW(hit) - n1 
PF = n1/NROW(hit)     # Percentage of Failures 

# Check if one of the series exhibits no-exceptions and then add one 
# violation in order to allow for estimation 
#if (n1 == 0 & PF ==0) { 
#   hit[NROW(hit)] = 1; 
#} 

# Kupiec (1995) Time Until First Failure 
if (n1 != 0 & PF !=0) { 
TUFF = head(which(hit==1),1) # Find the First Failure 
LRTUFF = -2*log((alpha*(1-alpha)^(TUFF-1))) + 2*log((1/TUFF)*(1-1/TUFF)^(TUFF-1)) # Log Likelihood of Time Until First Failure 
} else {  
TUFF = NA 
LRTUFF =NA} 

# Christoffersen (2003) Tests 
# Unconditional Coverage 
if (n1 != 0 & PF !=0) { 
        LRUC = -2*log(((alpha^n1)*((1-alpha)^n0))/((PF^n1)*(1-PF)^n0)) 
} else {LRUC = NA} 

# Independence Coverage 
if (n1 != 0 & PF !=0) { 
n00=n01=n10=n11=0 

for (i in 1:(NROW(fdata)-1)) { 
# Find 0 followed by 0 
        if (hit[i]==0 & hit[(i+1)] == 0) { n00 = n00 + 1} 
# Find 0 followed by 1 
        if (hit[i]==0 & hit[(i+1)] == 1) { n01 = n01 + 1} 
# Find 1 followed by 0 
        if (hit[i]==1 & hit[(i+0)] == 1) { n10 = n10 + 1} 
# Find 1 followed by 1 
        if (hit[i]==1 & hit[(i+1)] == 1) { n11 = n11 + 1} 
} 
rm(i) 

p01 = n01/(n00+n01) 
p00 = 1 - p01 
p11 = n11/(n10+n11) 
p2 = (n01+n11)/(n00+n01+n10+n11) 

# In case n11 = 0, then the test is estimated as ((1-p01)^n00)*(p01^n01) 
if (n11 == 0) { 
        LRIND = ((1-p01)^n00)*(p01^n01) 
                } else { 
    LRIND = -2*log((((1-p2)^(n00+n10))*(p2^(n01+n11)))/(((1-p01)^n00)*(p01^n01)*((1-p11)^n10)*(p11^n11))) 
} 
} else {LRIND = NA} 

# Conditional Coverage     
if (is.finite(LRUC)) {LRCC = LRUC + LRIND} else {LRCC = LRIND} 

# BASEL II Accord 
# According to Basel II, models are grouped into three categories: green, 
# yellow and red depending on the number of a% VaR violations.  
limits = cumsum(dbinom(1:50, size=NROW(fdata), prob=0.01)) 
green = length(which(limits<0.90)) 
yellow = green + length(which(limits>0.90 & limits<0.99)) 

if (n1 >= yellow) { 
        BASEL = -1 
        } else if (n1 <= yellow & n1 > green) { 
                BASEL = 0 
                } else { 
                BASEL = 1 
        } 

print(paste("Percentage of Failures: ", round(PF*100,2),"(%)",sep="")) 
print(paste("The Time of First Failure is: ", TUFF,sep=""))         
print(paste("Likelihood Ratio of Time of First Failure: ", round(LRTUFF,2),sep=""))         
print(paste("Likelihood Ratio of Unconditional Coverage: ", round(LRUC,2),sep=""))         
print(paste("Likelihood Ratio of Independence Coverage: ", round(LRIND,2),sep=""))         
print(paste("Likelihood Ratio of Conditional Coverage: ", round(LRCC,2),sep=""))         
print(paste("BASEL: ", BASEL,sep=""))         

# Results 
return = c(PF, TUFF, LRTUFF, LRUC, LRIND, LRCC, BASEL) 

} # End function 
