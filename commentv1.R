seir <- function(n=5.5e+6,ne=10,nt=150,gamma=1/3,delta=1/5) {
  # SEIR stochastic simulation model.
  # n = population size; ne = initially at E state; nt = number of days 
  # gamma = daily prob E -> I; delta = daily prob I -> R
  # expect the function outputs the simulation

  lambda <- 0.4/n
  x <- rep(0,n) # initialize to susceptible state
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) # individual infection rates and will use it to generate βi later

  # Now we want to find:
  # 1. the number of new infections each day
  # 2. the number of new infections among the 10% of the population with the lowest βi values
  # 3. the number of new infections in a random sample of 0.1% of the population.
  # So we firstly create list of index of the following features
  index1 <- order(beta)[1:(0.1*n)]# create a list of index of the number among the 10% of the population in increasing order using 'order' function  
  index2 <- sample(n,0.001*n,replace=FALSE)# create a list of index of the number in a random sample of 0.1% of the population
  x[1:ne] <- 1 # create some infectives

  # After creating the index lists, we now going to record the number of new infections of 3 different groups
  # (each day of whole population, 10% of the population with the lowest βi values, a random sample of 0.1% of the population)
  # by subtracting the number of elements in formal infection group from the number of elements in current infection group.
  # We need to calculate the number of elements of different groups,
  # computing the transmission process given in the context 
  # Here are the detailed and concrete steps
  S <- E_new <- E_low <- E_random <- S_low <- S_random <-rep(0,nt) ##set up storage for pop in each state 
  S[1] <- E_new[1] <- n-ne;
  S_low[1] <- sum(x[index1]==0); S_random[1] <- sum(x[index2]==0)
  E_low[1] <- sum(x[index1]==1); E_random[1] <- sum(x[index2]==1) ## initialize
  for (i in 2:nt) { ## loop over days
    u <- runif(n) ## uniform random deviates
    Iota<-sum(beta[x==2])
    x[x==2&u<delta] <- 3 ## I -> R with prob delta 
    x[x==1&u<gamma] <- 2 ## E -> I with prob gamma 
    x[x==0&u<lambda*beta*Iota] <- 1 ## S -> E with prob beta*I[i-1] 
    S[i] <- sum(x==0) #; E[i] <- sum(x==1); I[i] <- sum(x==2); R[i] <- sum(x==3)
    S_low[i]<- sum(x[index1]==0); S_random[i]<- sum(x[index2]==0)
    E_new[i]<- S[i-1]-S[i]; E_low[i]<- S_low[i-1]-S_low[i]; E_random[i]<- S_random[i-1]-S_random[i]
  }
  E_new[1] <- ne
  list(E_new=E_new,E_low=E_low,E_random=E_random) #E=E,I=I,R=R,
} ## seir

set.seed(3)
ep <- seir() ## run simulation 

n=5.5e+6
ep_new <- ep$E_new/n*10000
ep_low <- ep$E_low/n/0.1*10000
ep_random <- ep$E_random/n/0.001*10000

plot(ep_random,ylim=c(0,max(ep_new,ep_low,ep_random)),
     main="Daily Infection Trajectories",xlab="day",ylab="Incidence per 10000 per day", 
     type="l",col="dodgerblue",lwd=3) #
abline(h=max(ep_random), v=which(ep_random == max(ep_random)), lty=2, col="dodgerblue",lwd=2)
points(which(ep_random == max(ep_random)),max(ep_random), col="dodgerblue",cex=1.5,pch=16)
text(140,max(ep_random),  paste("peak at day",which(ep_random == max(ep_random))),col="brown",cex=1)

lines(ep_low,col="grey",lwd=4)
abline(h=max(ep_low), v=which(ep_low == max(ep_low)), lty=2, col= "grey",lwd=2)
points(which(ep_low == max(ep_low)),max(ep_low), col="grey",cex=1.5,pch=16)
text(140,max(ep_low), paste("peak at day",which(ep_low == max(ep_low))),col="brown",cex=1)

lines(ep_new, col="black",lwd=3)
abline(h=max(ep_new), v=which(ep_new == max(ep_new)), lty=2, col="black", lwd=2)
points(which(ep_new == max(ep_new)),max(ep_new), col="black",cex=1.5,pch=16)
text(140,max(ep_new), paste("peak at day",which(ep_new == max(ep_new))),col="brown",cex=1)
#zoom  

legend(list(x=2,y=140),legend=c("0.1% random","cautious 10%","whole population"),
       col=c("dodgerblue","grey","black"),lty=1,lwd=2,cex=0.9)

new <- low <- random <- list()
for (i in 1:10) {
  n=5.5e+6
  ep<-seir()
  new[[i]] <- ep$E_new
  low[[i]] <- ep$E_low
  random[[i]]<- ep$E_random
}


par(mfrow=c(3,1),mar=c(4,4,1,1),oma=c(0,0,0,0))
plot(new[[1]],xlab="day",ylab="Incidence per day",
     main = "Variability of 10 Repeated Simulations",
     type="l",col="black",lwd=2) #
legend("topleft","Whole Population",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")
max_new = mean(which(new[[1]]==max(new[[1]])))
for (i in 2:10) {
  lines(new[[i]],col="black",lwd=2)
  max_new  = max_new + mean(which(new[[i]]==max(new[[i]])))
}
abline(v=max_new/10, lty=2, col="red", lwd=2)
axis(side=1, at=max_new/10,labels=max_new/10,col ="red", col.axis="red")

plot(low[[1]],xlab="day",ylab="Incidence per day", type="l",col="grey",lwd=2) ## E black
legend("topleft","Cautious 10%",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")
max_low = mean(which(low[[1]]==max(low[[1]])))
for (i in 2:10) {
  lines(low[[i]],col="grey",lwd=2)
  max_low = max_low + mean(which(low[[i]]==max(low[[i]])))
}
abline(v=max_low/10, lty=2, col="red", lwd=2)
axis(side=1, at=max_low/10,labels=max_low/10,col="red",col.axis="red")

plot(random[[1]],xlab="day",ylab="Incidence per day", type="l",col="dodgerblue",lwd=2) #
legend("topleft","0.1% Random",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")
max_random = mean(which(random[[1]]==max(random[[1]])))
for (i in 2:10) {
  lines(random[[i]],col="dodgerblue",lwd=2)
  max_random = max_random + mean(which(random[[i]]==max(random[[i]])))
}
abline(v=max_random/10, lty=2, col="red", lwd=2)
axis(side=1, at=max_random/10,labels=max_random/10,col="red",col.axis="red")
