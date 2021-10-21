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
  # computing the transmission process given in the context and setting up storage for each state
  # Here are the detailed and concrete steps
  S <- E_new <- E_low <- E_random <- S_low <- S_random <-rep(0,nt) # set up storage for pop in each state 
  S[1] <- E_new[1] <- n-ne;# Using sum function to count the number of elements at 0-suscpetible state and initialize the number of people at susceptible state, removing the 10 E-state people
  S_low[1] <- sum(x[index1]==0)# Similarly, initialize the number of peiple at susceptible state among the 10% of the population with the lowest βi values
  S_random[1] <- sum(x[index2]==0)# initialize the number of people at susceptible state in a random sample of 0.1% of the population
  E_low[1] <- sum(x[index1]==1)# initialize the number of people at 1-exposed state among the 10% of the population with the lowest βi values
  E_random[1] <- sum(x[index2]==1) # initialize the number of people at exposed state in a random sample of 0.1% of the population
  for (i in 2:nt) { # loop over days
    u <- runif(n) # uniform random deviates
    Iota<-sum(beta[x==2])# create the set of all currently 2-infected people
    x[x==2&u<delta] <- 3 # I -> R with prob delta 
    x[x==1&u<gamma] <- 2 # E -> I with prob gamma 
    x[x==0&u<lambda*beta*Iota] <- 1 # S -> E with prob lambda*beta*Iota 
    S[i] <- sum(x==0)# assign the number of susceptible people each day of the whole population
    S_low[i]<- sum(x[index1]==0)# assign the number of susceptible people each day among the 10% of the population with the lowest βi values
    S_random[i]<- sum(x[index2]==0)# assign the number of susceptible people each day in a random sample of 0.1% of the population
    E_new[i]<- S[i-1]-S[i]# intepreting the number of new infections each day among the whole population by substracting the number of the day before from the number of current storage group
    E_low[i]<- S_low[i-1]-S_low[i]# intepreting the number of new infections each day among the 10% of the population with the lowest βi values by substracting the number of the day before from the number of current storage group
    E_random[i]<- S_random[i-1]-S_random[i]# Similar to above process
  }
  E_new[1] <- ne# initialize
  list(E_new=E_new,E_low=E_low,E_random=E_random) #E=E,I=I,R=R,
} ## seir

set.seed(3)
ep <- seir() # run simulation 

# Before constructing the plot, we standardize each trajectory by dividing by the size
n=5.5e+6
ep_new <- ep$E_new/n*10000 ## whole population
ep_low <- ep$E_low/n/0.1*10000 ## 10% of the population (with the lowest βi values)
ep_random <- ep$E_random/n/0.001*10000 ## 0.1% of the population

# Using plot function to plot the trajectory of new infections with peak in a random sample of 0.1% of the population
plot(ep_random,ylim=c(0,max(ep_new,ep_low,ep_random)),# set up the limit of the graph
     main="Daily Infection Trajectories",xlab="day",ylab="Incidence per 10000 per day", # write the title of the plot and lable axes
     type="l",col="dodgerblue",lwd=3) #using blue line to display the trajectory and assign the width
abline(h=max(ep_random), v=mean(which(ep_random == max(ep_random))), lty=2, col="dodgerblue",lwd=2)# add straight blue line with certain width to go through peak value and the corresponding line;using max and which to find the peak value and the corresponding day
points(mean(which(ep_random == max(ep_random))),max(ep_random), col="brown",cex=1.5,pch=16)# add a point of peak value of the trajectory
text(mean(which(ep_random == max(ep_random)))-2,max(ep_random)-6, mean(which(ep_random == max(ep_random))),col="brown",cex=1.2) # add explanation of "peak at day x", using which to find the certain day on which peak occurs

# Draw the trajectory of new infections with peak among the 10% of the population with the lowest βi values, the code is similar to above
lines(ep_low,col="grey",lwd=4)
abline(h=max(ep_low), v=mean(which(ep_low == max(ep_low))), lty=2, col= "grey",lwd=2)
points(mean(which(ep_low == max(ep_low))),max(ep_low), col="brown",cex=1.5,pch=16)
text(mean(which(ep_low == max(ep_low)))+2,max(ep_low)+6, mean(which(ep_low == max(ep_low))),col="brown",cex=1.2)

# Draw the trajectory of new infections with peak among the whole population, the code is similar to above
lines(ep_new, col="black",lwd=3)
abline(h=max(ep_new), v=mean(which(ep_new == max(ep_new))), lty=2, col="black", lwd=2)
points(mean(which(ep_new == max(ep_new))),max(ep_new), col="brown",cex=1.5,pch=16)
text(mean(which(ep_new == max(ep_new)))+2,max(ep_new)+6, mean(which(ep_new == max(ep_new))),col="brown",cex=1.2)
#zoom  

# Add legend indicating different groups to the plot, using list() to adjust the position of the legend to left of the plot
legend(list(x=2,y=140),legend=c("0.1% random","cautious 10%","whole population","day of peak"),
       col=c("dodgerblue","grey","black","brown"),lty=1,pch=c(NA,NA,NA,16),lwd=2,cex=0.9)

# Create 10 replicate simulations 
new <- low <- random <- list()
for (i in 1:10) {
  n=5.5e+6
  ep<-seir()
  new[[i]] <- ep$E_new
  low[[i]] <- ep$E_low
  random[[i]]<- ep$E_random
}

# Plot 10 repeated simulations of 3 groups in one plot and display the variability
# Firstly, draw trajectory for 'new infections each day among the whole population'
par(mfrow=c(3,1),mar=c(4,4,4,1),oma=c(0,0,0,0)) # Set margin and outer margin of the graph
plot(new[[1]],xlab="day",ylab="Incidence per day", # write title and label axex; Draw trajectory by lines
     main = "Variability of 10 Repeated Simulations",
     type="l",col="black",lwd=2) #
legend("topleft","Whole Population",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")# add legend
max_new = mean(which(new[[1]]==max(new[[1]])))# Find the peak and using mean function just in case there are more than one identical peak values
for (i in 2:10) {# loop
  lines(new[[i]],col="black",lwd=2)# draw new infections among whole population
  max_new  = max_new + mean(which(new[[i]]==max(new[[i]])))
}
abline(v=max_new/10, lty=2, col="red", lwd=2)# add line to indicate the peak
axis(side=1, at=max_new/10,labels=max_new/10,col ="red", col.axis="red")#make axes and lable them 

# Secondly, draw trajectory for 'the 10% of the population with the lowest βi values', similar process as above
plot(low[[1]],xlab="day",ylab="Incidence per day", type="l",col="grey",lwd=2) ## E black
legend("topleft","Cautious 10%",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")
max_low = mean(which(low[[1]]==max(low[[1]])))
for (i in 2:10) {
  lines(low[[i]],col="grey",lwd=2)
  max_low = max_low + mean(which(low[[i]]==max(low[[i]])))
}
abline(v=max_low/10, lty=2, col="red", lwd=2)
axis(side=1, at=max_low/10,labels=max_low/10,col="red",col.axis="red")

# Finally, draw trajectory for 'a random sample of 0.1% of the population', similar process as above
plot(random[[1]],xlab="day",ylab="Incidence per day", type="l",col="dodgerblue",lwd=2) #
legend("topleft","0.1% Random",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")
max_random = mean(which(random[[1]]==max(random[[1]])))
for (i in 2:10) {
  lines(random[[i]],col="dodgerblue",lwd=2)
  max_random = max_random + mean(which(random[[i]]==max(random[[i]])))
}
abline(v=max_random/10, lty=2, col="red", lwd=2)
axis(side=1, at=max_random/10,labels=max_random/10,col="red",col.axis="red")
