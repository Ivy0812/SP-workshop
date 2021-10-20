seir <- function(n=5.5e+6,ne=10,nt=150,gamma=1/3,delta=1/5) {
  ## SEIR stochastic simulation model.
  ## n = population size; ne = initially at E state; nt = number of days ## gamma = daily prob E -> I; delta = daily prob I -> R;
  lambda <- 0.4/n
  x <- rep(0,n) ## initialize to susceptible state
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) ## individual infection rates 
  #Now we want to find:
  #1. the number of new infections each day
  #2. the number of new infections among the 10% of the population with the lowest βi values
  #3. the number of new infections in a random sample of 0.1% of the population.
  #So we firstly create list of index with the following features
  index1 <- order(beta)[1:(0.1*n)]#create a list of index of the number of new infections among the 10% of the population in increasing order using 'order' function  
  index2 <- sample(n,0.001*n,replace=FALSE)
  x[1:ne] <- 1 ## create some infectives
  S <- E_new <- E_low <- E_random <- S_low <- S_random <-rep(0,nt) ## E <- I <- R <-  set up storage for pop in each state 
  S[1] <- E_new[1] <- n-ne;
  S_low[1] <- sum(x[index1]==0); S_random[1] <- sum(x[index2]==0)
  E_low[1] <- sum(x[index1]==1); E_random[1] <- sum(x[index2]==1) ## initialize
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    Iota<-sum(beta[x==1])
    x[x==2&u<delta] <- 3 ## I -> R with prob delta 
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma 
    x[x==0&u<lambda*beta*Iota] <- 1 ## S -> E with prob beta*I[i-1] 
    S[i] <- sum(x==0) #; E[i] <- sum(x==1); I[i] <- sum(x==2); R[i] <- sum(x==3)
    S_low[i]<- sum(x[index1]==0); S_random[i]<- sum(x[index2]==0)
    E_new[i]<- S[i-1]-S[i]; E_low[i]<- S_low[i-1]-S_low[i]; E_random[i]<- S_random[i-1]-S_random[i]
  }
  E_new[1] <- ne
  list(S=S,beta=beta,E_new=E_new,E_low=E_low,E_random=E_random) #E=E,I=I,R=R,
} ## seir
set.seed(3)
ep <- seir() ## run simulation 

n=5.5e+6
ep_new <- ep$E_new/n*10000
ep_low <- ep$E_low/n/0.1*10000
ep_random <- ep$E_random/n/0.001*10000

plot(ep_random,ylim=c(0,max(ep_new,ep_low,ep_random)),xlab="day",ylab="Incidence per 10000 day", type="l",col=4,lwd=2) ## E black
lines(ep_low,col="grey",lwd=4); lines(ep_new, col="black",lwd=3)

abline(h=max(ep_new), v=which(ep_new == max(ep_new)), lty=2, col="black", lwd=2)
abline(h=max(ep_low), v=which(ep_low == max(ep_low)), lty=2, col= "grey",lwd=2)
abline(h=max(ep_random), v=which(ep_random == max(ep_random)), lty=2, col=4,lwd=2)

points(which(ep_new == max(ep_new)),max(ep_new), col="black",cex=1.5,pch=16)
text(which(ep_new == max(ep_new))-12,max(ep_new)+3, "peak of whole population")

points(which(ep_low == max(ep_low)),max(ep_low), col="grey",cex=1.5,pch=16)
text(which(ep_low == max(ep_low))-12,max(ep_low)+3, "peak of cautious 10%")

points(which(ep_random == max(ep_random)),max(ep_random), col=4,cex=1.5,pch=16)
text(which(ep_random == max(ep_random))-12,max(ep_random)+3, "peak of 0.1% random sample")
#zoom

legend("top",legend=c("0.1% random","cautious 10%","whole population"),inset=c(0,-0.2),
       col=c(4,"grey","black"),lty=1,lwd=2,xpd = TRUE,horiz=T,cex=0.9)



plot(ep$E_new,ylim=c(0,100000),xlab="day",ylab="Incidence per 10000 day",type="l",col="black",lwd=1) ## E black

for (i in 1:9) {
  n=5.5e+6
  ep<-seir()
  ep_new <- ep$E_new
  lines(ep_new,col="black",lwd=1)
  abline(h=max(ep_new), v=which(ep_new == max(ep_new)), lty=2, col="black", lwd=1)
}

