install.packages("matrixStats")
library(matrixStats)

partitition=9901  #of quantiles here

teta=pi/partitition*0.5*(1:partitition)
N=7  #size of the first sample
M=13 
theta_known=exp(lnsigma)*sqrt(M/N)
dekning_faktisk
for