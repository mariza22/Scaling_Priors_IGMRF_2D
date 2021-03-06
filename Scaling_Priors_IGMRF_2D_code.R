##########################################################################
# This code provides all the steps that were used for producing the      #
# figures and the tables of the paper Scaling intrinsic Gaussian Markov  #
# random field priors in two dimensions. Describes the structure for     # 
# each type of IGMRFs in one and two dimensions.                         #
##########################################################################

library(ggplot2)

####  IGMRF/random walks ####

############ 1 dimension IGMRF #########################

## 1st order random walk, joint density (x1,x2,...,x100)
RW1 <- function(yr){
  R1                                                  <- matrix(0,yr,yr)
  R1[1,1:2] <- c(1,-1);R1[yr,(yr-1):yr]               <- c(-1,1)
  for (i in 2:(yr-1)) R1[i,(i-1):(i+1)]               <- c(-1,2,-1)
  eigenNoData                                         <- eigen(R1)
  eigenNoData$values[yr]                              <- Inf 
  EigenValGenInv                                      <- 1/eigenNoData$val
  EigenValGenInv[yr]                                  <- 0   
  SigmaGenInvNoTheta1                                 <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec))
  #plot(log(diag(SigmaGenInvNoTheta1)),xlab = "Node i", ylab="log-Marginal Standard Deviation") #with log 
  #plot(diag(SigmaGenInvNoTheta1),xlab = "Node i", ylab="Marginal Standard Deviation")
  return(sigma_ref_rw1 <- exp((1/yr)*sum(0.5*log(diag(SigmaGenInvNoTheta1)))) ) #page 5
}
## 2nd order random walk, joint density (x1,x2,...,x100)
RW2                                                   <- function(yr){
  R2                                                  <- matrix(0,yr,yr)
  R2[1,1:3] <- c(1,-2,1);R2[yr,(yr-2):yr]             <- c(1,-2,1)
  R2[2,1:4]                                           <- c(-2,5,-4,1);R2[yr-1,(yr-3):yr] <- c(1,-4,5,-2)
  for (i in 3:(yr-2))   R2[i,(i-2):(i+2)]             <- c(1,-4,6,-4,1)
  eigenNoData                                         <- eigen(R2)
  eigenNoData$values[(yr-1):yr]                       <- Inf 
  EigenValGenInv                                      <- 1/eigenNoData$val 
  EigenValGenInv[(yr-1):yr]                           <- 0  
  SigmaGenInvNoTheta2                                 <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec))
  #plot(log(diag(SigmaGenInvNoTheta2)),xlab = "Node i", ylab="log-Marginal Standard Deviation")
  #plot(diag(SigmaGenInvNoTheta2),xlab = "Node i", ylab="Marginal Standard Deviation")
  return(sigma_ref_rw2 <- exp((1/yr)*sum(0.5*log(diag(SigmaGenInvNoTheta2)))) )
}
#### 2 dimensions IGMRF ####

## 2D second order random walk, (bound 1, Yue & Speckman )
RW2D                                                  <- function(yr){
  mtx3 <- diag(1,yr) ;mtx1 <- mtx2 <- mtx4 <- mtx5<- mtx6 <- matrix(0,yr,yr)
  # 2nd order random walk
  mtx1[1,1:3] <- mtx1[yr,yr:(yr-2)]                   <- c(6,-5,1)  
  mtx1[2,1:4] <- mtx1[(yr-1),yr:(yr-3)]               <- c(-5,12,-6,1) 
  for (i in 3:(yr-2)) mtx1[i,(i-2):(i+2)]             <- c(1,-6,12,-6,1) 
  #1st order random walk
  mtx2[1,1:2] <-  mtx2[yr,yr:(yr-1)]                  <- c(-5,2)
  for (i in 2:(yr-1)) mtx2[i,(i-1):(i+1)]             <- c(2,-7,2)
  # 2nd order random walk
  mtx4[1,1:3] <- mtx4[yr,yr:(yr-2)]                   <- c(12,-7,1)
  mtx4[2,1:4] <- mtx4[(yr-1),yr:(yr-3)]               <- c(-7,20,-8,1)
  for (i in 3:(yr-2)) mtx4[i,(i-2):(i+2)]             <- c(1,-8,20,-8,1) 
  #1st order random walk
  mtx5[1,1:2] <-  mtx5[yr,yr:(yr-1)]                  <- c(-6,2)
  for (i in 2:(yr-1)) mtx5[i,(i-1):(i+1)]             <- c(2,-8,2)
  R2D                                                 <- matrix(0,yr^2,yr^2)
  R2D[1:yr,1:(3*yr)]                                  <- cbind(mtx1,mtx2,mtx3)
  R2D[((yr-1)*yr+1):(yr^2),((yr-3)*yr+1):(yr^2)]      <- cbind(mtx3,mtx2,mtx1)
  R2D[(yr+1):(2*yr),1:(4*yr)]                         <- cbind(mtx2,mtx4,mtx5,mtx3)
  R2D[((yr-2)*yr+1):((yr-1)*yr),((yr-4)*yr+1):(yr^2)] <- cbind(mtx3,mtx5,mtx4,mtx2)
  for ( i in 3:(yr-2)) R2D[((i-1)*yr+1):(i*yr),((i-3)*yr+1):((i+2)*yr)] <- cbind(mtx3,mtx5,mtx4,mtx5,mtx3)
  #print(qr(R2D)$rank)
  eigenNoData                                         <- eigen(R2D)
  eigenNoData$values[(yr^2-2):(yr^2)]                 <- Inf 
  EigenValGenInv                                      <- 1/eigenNoData$val
  EigenValGenInv[(yr^2-2):(yr^2)]                     <- 0  
  SigmaGenInvNoTheta3                                 <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec))
  #plot(log(diag(SigmaGenInvNoTheta3)),xlab = "Node i", ylab="Marginal Standard Deviation")
  return(sigma_ref_rw2D <- exp((1/(yr^2))*sum(0.5*log(diag(SigmaGenInvNoTheta3)))) )
}

#Figure 1

RW1(50)
RW2(50)
RW2D(50) #it takes time


## 2D Rue&Held, 1st suggestion (Torus 1)
RW2D_1                                                <- function(yr){
  b1                  <- diag(-4, yr);b2              <- diag(1,yr)
  b3                  <- matrix(0,yr,yr);D            <- matrix(0,yr^2,yr^2)
  b1[abs(row(b1) - col(b1)) == 1]                     <- 1
  b1[1,yr]            <- b1[yr,1]                     <- 1
  for ( i in 1:yr) D[((i-1)*yr+1): (i*yr),((i-1)*yr+1):(i*yr)] <- b1
  for ( i in 1:(yr-1)){
    D[((i-1)*yr+1):(i*yr),(i*yr+1):(i*yr+yr)]         <- b2
    D[(i*yr+1):(i*yr+yr),((i-1)*yr+1):(i*yr)]         <- b2
  }
  D[1:yr,(yr*(yr-1)+1):(yr^2)]<-D[(yr*(yr-1)+1):(yr^2),1:yr]<- b2
  D1                                                  <- t(D)%*%D
  #print(qr(D1)$rank)  # since we are in a torus the rank defficiency is 1 not 3
  # qr(D)$rank        #n^2-1
  eigenNoData                                         <- eigen(D1)
  eigenNoData$values[(yr^2-2):(yr^2)]                 <- Inf 
  EigenValGenInv                                      <- 1/eigenNoData$val
  EigenValGenInv[(yr^2-2):(yr^2)]                     <- 0  
  SigmaGenInvNoTheta31                                <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec))
  plot(log(diag(SigmaGenInvNoTheta31)))
  return( exp((1/(yr^2))*sum(0.5*log(diag(SigmaGenInvNoTheta31)))) )
  #print(D1)
}
##  2D Rue&Held, 2nd suggestion  (Torus 2)
RW2D_2                                                <- function(yr){ 
  b1 <- diag(-20,yr);b2                               <- diag(4,yr)
  b3 <- matrix(0,yr,yr);D                             <- matrix(0,yr^2,yr^2)
  b1[abs(row(b1) - col(b1)) == 1]                     <- 4
  b2[abs(row(b2) - col(b2)) == 1]                     <- 1
  b1[1,yr]               <- b1[yr,1]                  <- 4
  b2[1,yr]               <- b2[yr,1]                  <- 1
  for ( i in 1:yr) D[((i-1)*yr+1):(i*yr),((i-1)*yr+1):(i*yr)] <- b1
  for ( i in 1:(yr-1)){
    D[((i-1)*yr+1): (i*yr), (i*yr+1):(i*yr+yr)]       <- b2
    D[(i*yr+1):(i*yr+yr), ((i-1)*yr+1):(i*yr)]        <- b2
  }
  D[1:yr,(yr*(yr-1)+1):(yr^2)] <- D[(yr*(yr-1)+1):(yr^2),1:yr]<- b2
  D2                                                  <- t(D)%*%D
  #qr(D)$rank  # since we are in a torus the rank defficiency is 1 not 3!
  #print(qr(D2)$rank)
  eigenNoData                                         <- eigen(D2)
  eigenNoData$values[(yr^2-2):(yr^2)]                 <- Inf 
  EigenValGenInv                                      <- 1/eigenNoData$val
  EigenValGenInv[(yr^2-2):(yr^2)]                     <- 0  
  SigmaGenInvNoTheta32                                <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec))
  plot(log(diag(SigmaGenInvNoTheta32)))
  return( exp((1/yr^2)*sum(0.5*log(diag(SigmaGenInvNoTheta32)))) )
}
## 2D second order random walk, (bound 2, Terzopoulos )
RW2D_terz                                             <- function(yr){
  mtx3 <- diag(1,yr) ;mtx1 <- mtx2 <- mtx4 <- mtx5<- mtx6 <- matrix(0,yr,yr)
  # 2nd order random walk
  mtx1[1,1:3] <- mtx1[yr,yr:(yr-2)]                   <- c(4,-4,1)  
  mtx1[2,1:4] <- mtx1[(yr-1),yr:(yr-3)]               <- c(-4,10,-6,1) 
  for (i in 3:(yr-2)) mtx1[i,(i-2):(i+2)]             <- c(1,-6,11,-6,1) 
  #1st order random walk
  mtx2[1,1:2] <-  mtx2[yr,yr:(yr-1)]                  <- c(-4,2)
  for (i in 2:(yr-1)) mtx2[i,(i-1):(i+1)]             <- c(2,-6,2)
  # 2nd order random walk
  mtx4[1,1:3] <- mtx4[yr,yr:(yr-2)]                   <- c(10,-6,1)
  mtx4[2,1:4] <- mtx4[(yr-1),yr:(yr-3)]               <- c(-6,18,-8,1)
  for (i in 3:(yr-2)) mtx4[i,(i-2):(i+2)]             <- c(1,-8,19,-8,1) 
  #1st order random walk
  mtx5[1,1:2] <-  mtx5[yr,yr:(yr-1)]                  <- c(-6,2)
  for (i in 2:(yr-1)) mtx5[i,(i-1):(i+1)]             <- c(2,-8,2)
  # 2nd order random walk
  mtx6[1,1:3] <- mtx6[yr,yr:(yr-2)]                   <- c(11,-6,1)
  mtx6[2,1:4] <- mtx6[(yr-1),yr:(yr-3)]               <- c(-6,19,-8,1)
  for (i in 3:(yr-2)) mtx6[i,(i-2):(i+2)]             <- c(1,-8,20,-8,1) 
  
  R2D                                                 <- matrix(0,yr^2,yr^2)
  R2D[1:yr,1:(3*yr)]                                  <- cbind(mtx1,mtx2,mtx3)
  R2D[((yr-1)*yr+1):(yr^2),((yr-3)*yr+1):(yr^2)]      <- cbind(mtx3,mtx2,mtx1)
  R2D[(yr+1):(2*yr),1:(4*yr)]                         <- cbind(mtx2,mtx4,mtx5,mtx3)
  R2D[((yr-2)*yr+1):((yr-1)*yr),((yr-4)*yr+1):(yr^2)] <- cbind(mtx3,mtx5,mtx4,mtx2)
  for ( i in 3:(yr-2)) R2D[((i-1)*yr+1):(i*yr),((i-3)*yr+1):((i+2)*yr)] <- cbind(mtx3,mtx5,mtx6,mtx5,mtx3)
  
  eigenNoData                                         <- eigen(R2D)
  eigenNoData$values[(yr^2-2):(yr^2)]                 <- Inf 
  EigenValGenInv                                      <- 1/eigenNoData$val##!!!!! 
  EigenValGenInv[(yr^2-2):(yr^2)]                     <- 0  
  SigmaGenInvNoTheta3                                 <- eigenNoData$vec %*% (EigenValGenInv * t(eigenNoData$vec))
  plot(log(diag(SigmaGenInvNoTheta3)),xlab = "Node i", ylab="Marginal Standard Deviation")
  return(sigma_ref_rw2D <- exp((1/(yr^2))*sum(0.5*log(diag(SigmaGenInvNoTheta3)))) )
}
#### Application of the  Sections 3 and 4 #### 
## Section 3 ##

# FIGURE 2
s_RW1                                                 <- s_RW2 <- s_RW2D <- s_RW2D_terz <- s_RW2D_1<-s_RW2D_2 <-rep(0,12)
nodes   <- seq(5,60,5);k                              <- 1
for (i in nodes){
  s_RW1[k]                                            <- RW1(yr=i)  
  s_RW2[k]                                            <- RW2(yr=i) 
  s_RW2D[k]                                           <- RW2D(yr=i)
  s_RW2D_terz[k]                                      <- RW2D_terz(yr=i)
  s_RW2D_1[k]                                         <- RW2D_1(yr=i) 
  s_RW2D_2[k]                                         <- RW2D_2(yr=i) 
  k                                                   <- k+1
}
plot(nodes,s_RW1,type='l',ylim=c(0,20),lty=2,xlab = "Number of nodes",ylab="Standard Deviation")
lines(nodes,s_RW2,col=2,lty=3)
lines(nodes,s_RW2D,col=3,lty=4)
lines(nodes,s_RW2D_terz,col=4,lty=4)
legend(5, 20, c("RW1","RW2","RW2D","RW2D2"), col = 1:4, lty = 2:4,ncol = 3)

df                                                    <- data.frame(Nodes = seq(5,60,5),
                                                                    Standard_Deviation = c(s_RW1,s_RW2,s_RW2D),#,s_RW2D_terz
                                                                    Type = rep(c("s_RW1","s_RW2","s_RW2D"),each=12))#,"s_RW2D_terz"

gg                                                    <- ggplot(data=df, aes(x=Nodes,y=Standard_Deviation,group=Type))+
  geom_line(aes(color=Type))+
  geom_point(aes(color=Type))+
  labs(title="",x="Nodes", y = "Standard Deviation")
gg + theme(text = element_text(size = 15))    

# TABLE 1
for ( i in c(11,20,40)){#,100
  print(RW2D_1(i) )
  print(RW2D_2(i) )
  print( RW2D(i)  )
  print(RW2D_terz(i) )
}
## Section 4 ##

yr                                                    <- 40
s_RW2                                                 <- RW2(yr=yr)
s_RW2D                                                <- RW2D(yr=yr) 
mu                                                    <- 7         
b                                                     <- 2
alpha                                                 <- 0.001

#upper limits of the table 10. The U formula is used

# Normal distribution
U_rw2                                                 <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=mu,sd=1))
U_rw2D                                                <- sqrt( b*(s_RW2D^2)/qnorm(alpha,mean=mu,sd=1))
Us1                                                   <- c(U_rw2,U_rw2D)
summary(Us1)                     #Based on median we choose u=0.5
u1                                                    <- round(summary(Us1)[[3]],3)
#I will use the formula of U in respect of beta taking the median of the summary(u)
b_rw2                                                 <- (u1^2)*(qnorm(alpha, mean=mu, sd=1))/(s_RW2^2)
b_rw2D                                                <- (u1^2)*(qnorm(alpha, mean=mu, sd=1))/(s_RW2D^2)
print(c(b_rw2,b_rw2D))

#TABLE 2
s_RW11                                                <- RW1(yr=11)
s_RW12                                                <- RW2(yr=11)
s_RW12D                                               <- RW2D(yr=11) 
s_RW21                                                <- RW1(yr=20)
s_RW22                                                <- RW2(yr=20)
s_RW22D                                               <- RW2D(yr=20) 

#TABLE 3
for ( i in c(1,2,3)){
  for ( j in c(5,11,20,40)){
    print(paste("b is: ",i))  
    print(paste("years are: ",j)) 
    s_RW1                                              <- RW1(yr=j) 
    s_RW2                                              <- RW2(yr=j)
    s_RW2D                                             <- RW2D(yr=j) 
    mu                                                 <- 7         
    b                                                  <- i
    alpha                                              <- 0.001
    U_rw1                                              <- sqrt( b*(s_RW1^2)/qnorm(alpha,mean=mu,sd=1))
    U_rw2                                              <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=mu,sd=1))
    U_rw2D                                             <- sqrt( b*(s_RW2D^2)/qnorm(alpha,mean=mu,sd=1))
    Us1                                                <- c(U_rw1,U_rw2,U_rw2D)
    summary(Us1) #Based on median we choose u=0.5
    u                                                  <- round(summary(Us1)[[3]],3)
    b_rw1                                              <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW1^2)
    b_rw2                                              <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2^2)
    b_rw2D                                             <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2D^2)
    print(c(b_rw1,b_rw2,b_rw2D))
  }
}
#TABLE 4 

k                <- 1;       m                         <- c(7,8,10,12)
for ( i in c(0.9,1.2,1.59,3.55)){
  j<- 11
  print(paste("b is: ",i))  
  s_RW2                                                <- RW2(yr=j)
  s_RW2D                                               <- RW2D(yr=j) 
  mu                                                   <- m[k]       
  b                                                    <- i
  alpha                                                <- 0.001
  U2_rw2                                               <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=mu,sd=1))#abs()
  U2_rw2D                                              <- sqrt( b*(s_RW2D^2)/qnorm(alpha,mean=mu,sd=1))#abs()
  Us1                                                  <- c(U2_rw2,U2_rw2D)
  summary(Us1) 
  u                                                    <- round(summary(Us1)[[3]],3)
  b2_rw2                                               <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2^2)
  b2_rw2D                                              <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2D^2)
  print(c(b2_rw2,b2_rw2D))
  k<-k+1
}

###### EXTRA: what if I change the U value and the alpha parameter ? in the table 4. ####

#alteration of the 4 table 
# change of alpha, a=0.01

for ( i in c(0.9,1.2,1.59,3.55)){
  j<- 11
  print(paste("b is: ",i))  
  s_RW2         <- RW2(yr=j)
  s_RW2D        <- RW2D(yr=j) 
  mu         <- 7         
  b             <- i
  alpha         <- 0.01
  U2_rw2        <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=mu,sd=1))
  U2_rw2D       <- sqrt( b*(s_RW2D^2)/qnorm(alpha,mean=mu,sd=1))
  
  Us1           <- c(U2_rw2,U2_rw2D)
  summary(Us1) #Based on median we choose u=0.5
  u             <- round(summary(Us1)[[3]],3)
  b2_rw2        <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2^2)
  b2_rw2D       <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2D^2)
  print(c(b2_rw2,b2_rw2D))
}
#  0.5334365 1.8354918
#  0.7098331 2.4424515
#  0.9410647 3.2380921
#  2.101108 7.229664

# we increase the probability in the tails and the outcome seems more less the same as the table 4.

#alteration of the 4 table 
# change of U: 3rd quantile, 5th position
for ( i in c(0.9,1.2,1.59,3.55)){
  j<- 11
  print(paste("b is: ",i))  
  s_RW2         <- RW2(yr=j)
  s_RW2D        <- RW2D(yr=j) 
  mu         <- 7         
  b             <- i
  alpha         <- 0.001
  U2_rw2        <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=mu,sd=1))
  U2_rw2D       <- sqrt( b*(s_RW2D^2)/qnorm(alpha,mean=mu,sd=1))
  
  Us1           <- c(U2_rw2,U2_rw2D)
  summary(Us1) #Based on median we choose u=0.5
  u             <- round(summary(Us1)[[5]],3) # 3rd quantile, 5th position
  b2_rw2        <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2^2)
  b2_rw2D       <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2D^2)
  print(c(b2_rw2,b2_rw2D))
}

#alteration of the 4 table 
# change of U #1st quantile, 2nd position
for ( i in c(0.9,1.2,1.59,3.55)){
  j<- 11
  print(paste("b is: ",i))  
  s_RW2         <- RW2(yr=j)
  s_RW2D        <- RW2D(yr=j) 
  mu         <- 7         
  b             <- i
  alpha         <- 0.001
  U2_rw2        <- sqrt( b*(s_RW2^2)/qnorm(alpha,mean=mu,sd=1))
  U2_rw2D       <- sqrt( b*(s_RW2D^2)/qnorm(alpha,mean=mu,sd=1))
  
  Us1           <- c(U2_rw2,U2_rw2D)
  summary(Us1) #Based on median we choose u=0.5
  u             <- round(summary(Us1)[[2]],3) # 3rd quantile, 5th position
  b2_rw2        <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2^2)
  b2_rw2D       <- (u^2)*(qnorm(alpha,mean=mu, sd=1))/(s_RW2D^2)
  print(c(b2_rw2,b2_rw2D))
}

#the values for the variance for the 2d model can change a lot with the U upper level. 
#so if you want smaller values than those in the 4th table use the U of the 1st quantile.
