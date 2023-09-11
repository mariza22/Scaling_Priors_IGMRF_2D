#Paper analysis: Scaling intrinsic Gaussian Markov random field priors in spatial modelling
###### 1st order random walk, joint density (x1,x2,...x100)
RW1 <- function(yr){                                                            # number of nodes (or years) 
  R1                             <- matrix(0,yr,yr)                             # defining the precision matrix 
  R1[1,1:2]                      <- c(1,-1);R1[yr,(yr-1):yr]<-c(-1,1)           # first and last row 
  for (i in 2:(yr-1)) R1[i,(i-1):(i+1)] <- c(-1,2,-1)                           # central rows 
  eigenNoData                    <- eigen(R1)                                   # eigenvalues and eigenvector 
  #eigenNoData$values[yr]         <- Inf 
  EigenValGenInv                 <- 1/eigenNoData$val                           # the inverse of the eigenvalues
  EigenValGenInv[yr]             <- 0                                           # setting the inverse eigenvalues =0 too
  SigmaGenInvNoTheta1            <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec)) # calculating the Sigma* matrix
  print(sigma_ref_rw1            <- exp((1/yr)*sum(0.5*log(diag(SigmaGenInvNoTheta1)))) )     # calculating the Sigma reference
  #plot(log(diag(SigmaGenInvNoTheta1)),xlab = "Node i", ylab="Marginal Standard Deviation")
}
RW1(yr=10)
###### 2nd order random walk, joint density (x1,x2,...x100)
RW2<- function(yr){                                                             # number of nodes (or years)
  R2                             <- matrix(0,yr,yr)                             # defining the precision matrix
  R2[1,1:3]                      <- c(1,-2,1);R2[yr,(yr-2):yr] <- c(1,-2,1)     # 1st and last row
  R2[2,1:4]                      <- c(-2,5,-4,1);R2[yr-1,(yr-3):yr] <- c(1,-4,5,-2) # 2nd and 2nd last row
  for (i in 3:(yr-2))   R2[i,(i-2):(i+2)] <- c(1,-4,6,-4,1)                     # central rows
  eigenNoData                    <- eigen(R2)                                   # eigenvalues and eigenvector
  #eigenNoData$values[(yr-1):yr]  <- Inf                                        
  EigenValGenInv                 <- 1/eigenNoData$val                           # the inverse of the eigenvalues 
  EigenValGenInv[(yr-1):yr]      <- 0                                           # setting the inverse eigenvalues =0 too
  SigmaGenInvNoTheta2            <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec)) # calculating the Sigma* matrix
  print(sigma_ref_rw2            <- exp((1/yr)*sum(0.5*log(diag(SigmaGenInvNoTheta2)))) )     # calculating the Sigma reference
  #plot(log(diag(SigmaGenInvNoTheta2)),xlab = "Node i", ylab="Marginal Standard Deviation")
}
RW2(yr=10)

#Reproducing Matrix 10: The sigma references can occur when 
#yr <- 496 ; yr<- 16 for the rw1 and rw2( code line: 2-26) for wbc and tpi variables respectively. 
RW1(yr=496); RW2(yr=496) #for wbc
RW1(yr=16);  RW2(yr=16) #for tpi


s_besag                          <- 0.64

# Sigma reference for a different number of nodes (years) of 1st and 2nd order random walk

s_wbc_rw1                        <- RW1(yr=496)                    
s_wbc_rw2                        <- RW2(yr=496)
s_tpi_rw1                        <- RW1(yr=16)
s_tpi_rw2                        <- RW2(yr=16)

#Setting initial values
shape                            <- 1
b                                <- 5*10^(-5)  #initial scalar, variance parameter
alpha                            <- 0.001

# Upper limits of the table 10. The U formula is used
# We calculate the upper limits for the 1st, 2nd order random walk

U_besag                          <- sqrt( b*(s_besag^2)/(qgamma(alpha,shape=shape,rate=1)))
U_wbc_rw1                        <- sqrt( b*(s_wbc_rw1^2)/qgamma(alpha,shape=shape,rate=1))
U_wbc_rw2                        <- sqrt( b*(s_wbc_rw2^2)/qgamma(alpha,shape=shape,rate=1))
U_tpi_rw1                        <- sqrt( b*(s_tpi_rw1^2)/qgamma(alpha,shape=shape,rate=1))
U_tpi_rw2                        <- sqrt( b*(s_tpi_rw2^2)/qgamma(alpha,shape=shape,rate=1))
us                               <- c(U_besag,U_wbc_rw1,U_wbc_rw2,U_tpi_rw1,U_tpi_rw2)
summary(us)                             #Based on median we choose u=0.5
u                                <- 0.1 #summary(us)[3]

#I will use the formula of U in respect of beta

(b_besag                          <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_besag^2))
(b_wbc_rw1                        <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_wbc_rw1^2))
(b_wbc_rw2                        <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_wbc_rw2^2))
(b_tpi_rw1                        <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_tpi_rw1^2))
(b_tpi_rw2                        <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_tpi_rw2^2))

################# changing distribution  ###################
##### using number of nodes 16 and 60. 

s_RW1.1                          <- RW1(yr=16) 
s_RW1.2                          <- RW1(yr=60) 
s_RW2.1                          <- RW2(yr=16) 
s_RW2.2                          <- RW2(yr=60)

shape                            <- 1               # shape is used for the gamma distribution of the F-1(.)
b                                <- 5*10^(-5)       # different value for the initial standard deviation
alpha                            <- 0.001

#upper limits of the table 10. The U formula is used
###Gamma distribution
U_rw1.1                          <- sqrt( b*(s_RW1.1^2)/qgamma(alpha,shape=shape,rate=1))
U_rw1.2                          <- sqrt( b*(s_RW1.2^2)/qgamma(alpha,shape=shape,rate=1))
U_rw2.1                          <- sqrt( b*(s_RW2.1^2)/qgamma(alpha,shape=shape,rate=1))
U_rw2.2                          <- sqrt( b*(s_RW2.2^2)/qgamma(alpha,shape=shape,rate=1))
Us                               <- c(U_rw1.1,U_rw1.2,U_rw2.1,U_rw2.2)
summary(Us) #Based on median we choose u=0.5
u                                <- round(summary(Us)[[3]],3)     # choice of Upper level
### Normal distribution    ####
# Instead of using a gamma distribution we will proceed to a truncated normal disribution

alpha1                           <- 0.001                         # same as alpha
U2_rw1.1                         <- sqrt( b*(s_RW1.1^2)/abs(qnorm(alpha1,mean=shape,sd=1)))
U2_rw1.2                         <- sqrt( b*(s_RW1.2^2)/abs(qnorm(alpha1,mean=shape,sd=1)))
U2_rw2.1                         <- sqrt( b*(s_RW2.1^2)/abs(qnorm(alpha1,mean=shape,sd=1)))
U2_rw2.2                         <- sqrt( b*(s_RW2.2^2)/abs(qnorm(alpha1,mean=shape,sd=1)))
Us1                              <- c(U2_rw1.1,U2_rw1.2,U2_rw2.1,U2_rw2.2)
summary(Us1) #Based on median we choose u=0.5
u1                               <- round(summary(Us1)[[3]],3) # choice of Upper level using the median

#I will use the formula of U in respect of beta
#taking the median of the summary(u)/summary(u1)

#gamma distribution
b_rw1.1                          <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_RW1.1^2)
b_rw1.2                          <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_RW1.2^2)
b_rw2.1                          <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_RW2.1^2)
b_rw2.2                          <- (u^2)*(qgamma(alpha,shape=shape, rate=1))/(s_RW2.2^2)
print(c(b_rw1.1,b_rw1.2,b_rw2.1,b_rw2.2))

#normal distribution
b2_rw1.1                         <- (u1^2)*(abs(qnorm(alpha1,mean=shape, sd=1))/(s_RW1.1^2))
b2_rw1.2                         <- (u1^2)*(abs(qnorm(alpha1,mean=shape, sd=1))/(s_RW1.2^2))
b2_rw2.1                         <- (u1^2)*(abs(qnorm(alpha1,mean=shape, sd=1))/(s_RW2.1^2))
b2_rw2.2                         <- (u1^2)*(abs(qnorm(alpha1,mean=shape, sd=1))/(s_RW2.2^2))
print(c(b2_rw1.1,b2_rw1.2,b2_rw2.1,b2_rw2.2))

##### different prior, application to our problem  ###### 

#defining the prior 
nds                              <- seq(5,30,5)               # the number of nodes vary from 5 to 30
shape                            <- 8                         # the prior: Normal(7,1.25)
b                                <- 5                         # initial value of standard deviation
alpha                            <- .01
b_rw1_store                      <- rep(0,length(nds) )       # the new std for random walk order 1 depending on the number of nodes
b_rw2_store                      <- rep(0,length(nds) )       # the new std for random walk order 2 depending on the number of nodes
for (i in 1:length(nds)) {                                    # position of number of nodes
    s_RW1                        <- RW1(yr=nds[i])            # sigma reference for RW1 
    s_RW2                        <- RW2(yr=nds[i])            # sigma reference for RW2
    # The U formula is used
    (U_rw1                       <- sqrt( b*(s_RW1^2)/qnorm(alpha,mean=shape,sd=1)) )# calculating upper levels
    (U_rw2                       <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=shape,sd=1)) )
    Us                           <- c(U_rw2,U_rw1) 
    u                            <- round(summary(Us)[[3]],3)#upper[j]# choosing the median of the upper levels 
    b_rw1_store[i]               <- (u^2)*(qnorm(alpha,mean=shape, sd=1))/(s_RW1^2) #the suggested varia 
    b_rw2_store[i]               <- (u^2)*(qnorm(alpha,mean=shape, sd=1))/(s_RW2^2)
  }

#different values for b

b_options                        <- seq(0.5,4,0.05)

for (j in b_options ){
  print(paste("b:",j))
  i                              <- 11#c(5,10,20,40)
  print(paste("years are:",i))
  s_RW1                          <- RW1 (yr=i)  # 5 years 
  s_RW2                          <- RW2 (yr=i) 
  shape                          <- 7           #the prior: Normal(7,1.25)
  b                              <- j 
  alpha                          <- .001
  # The U formula is used
  U_rw1                          <- sqrt( b*(s_RW1^2)/qnorm(alpha,mean=shape,sd=1))
  U_rw2                          <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=shape,sd=1))
  Us                             <- c(U_rw1, U_rw2)
  summary(Us) 
  u                              <- round(summary(Us)[[3]],3)
  #I will use the formula of U in respect of beta
  b_rw1                          <- (u^2)*(qnorm(alpha,mean=shape, sd=1))/(s_RW1^2)
  print(paste("b_rw1:", b_rw1) )
  b_rw2                          <- (u^2)*(qnorm(alpha,mean=shape, sd=1))/(s_RW2^2)
  print(paste("b_rw2:", b_rw2) )
}
#####
#known the scaling of the 1st order random walk we can calculate easily the scaling for the 2nd order random walk 
b_rw1                            <- 3
(b_rw2                           <- b_rw1*(s_RW1^2/s_RW2^2))

#known the scaling of the 2nd order random walk we can calculate easily the scaling for the 1st order random walk 
b_rw2                            <- 3
(b_rw1                           <- b_rw2*(s_RW2^2/s_RW1^2)) 
