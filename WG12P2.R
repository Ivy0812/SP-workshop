# This work is completed by group 12：Wenyi Yu s2161093, Minke Pan s2160782 and Yanren Mao s2207399
# address of github repo: https://github.com/Ivy0812/SP-workshop.git

# Overview
# In order to solve these 2 questions:
# 1.What do the daily infection trajectories look like for the population with the lowest individual transmission probabilities?
# 2.How do they compare to the infection trajectories for the whole population?

# Firstly we create a SEIR (Susceptible-Exposed-Infective-Recovered) models.
# Almost everyone starts out susceptible, and then move to the exposed class 
# according to their daily probability of being infected by someone in the infected class. 
# In particular suppose that the ith person in the whole population has a relative contact rate with other people of βi.
# Each day an uninfected person, j, has a probability λ*βj*sum(βi) of being infected, and entering the exposed (E) state.
# People move from the exposed to infected class with a constant probability, gamma , each day, 
# and move from the infected to the R class with a different constant probability, delta.
# Suppose that we want to simulate from this stochastic SEIR model, maintaining a record of the new infections in the different groups over time. 
# We produce a function that takes parameter values and control constants as inputs, and outputs vectors of the new infections of the whole, 10% and 0.1% populations over time.
# We create a single vector maintaining an integer state indicating each individual’s current status 0 for S, 1 for E, 2 for I and 3 for R.
# The individual states are then updated daily according to the model, and the numbers in each state summed.

# Secondly, We standardize each trajectory by dividing their population size to plot them on the same plot.
# Then we produce a plot showing how the daily infection trajectories compare between 
# the whole population, the ‘cautious 10%’ and the 0.1% random sample 
# and also provide labels giving the day on which each trajectory peaks.

# Finally we write code to visualize the variability in the results from running 10 replicate simulations，
# and record the 10 daily infection trajectories of each group by plotting 3 plots respectively.
# The last plot shows the difference of the day on which the peak value occurs between whole population, the ‘cautious 10%’ and the 0.1% random sample. 


# Step 1: Simulation
seir <- function(n=5.5e+6,ne=10,nt=150,gamma=1/3,delta=1/5) {
  # SEIR stochastic simulation model.
  # n = population size; ne = initially at exposed state; nt = number of days 
  # gamma = daily prob E -> I; delta = daily prob I -> R
  # expect the function outputs the simulation

  lambda <- 0.4/n
  x <- rep(0,n) # initialize to susceptible state
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) # individual infection rates and will use it to generate βi later

  # Now we want to find:
  # 1. the number of new infections each day among the whole population
  # 2. the number of new infections among the 10% of the population with the lowest βi values
  # 3. the number of new infections in a random sample of 0.1% of the population.

  index1 <- order(beta)[1:(0.1*n)]# fix indices of the lowest 10% βi values of the population
  index2 <- sample(n,0.001*n,replace=FALSE)# randomly sample indices of 0.1% of the population
  

  # After creating the index lists, we now going to track the number of new infections of 3 different groups:
  # 1.each day of whole population
  # 2.10% of the population with the lowest βi values
  # 3.a random sample of 0.1% of the population
  # We calculate the new infections of each group by computing the decreases of the susceptible state between 2 consecutive days
  # Here are the simulation of the transmission process
  
  S <- E_new <- E_low <- E_random <- S_low <- S_random <-rep(0,nt) # set up storage for pop in each state 
  x[1:ne] <- 1; E_new[1] <- ne  # setting 10 chosen people in exposed state
  S[1] <- n-ne# removing 10 E-state people from S state
  
  # assign the number of people at S state and E state among the 10% of the population
  S_low[1] <- sum(x[index1]==0) ; E_low[1] <- sum(x[index1]==1)
  # assign the number of people at S state and E state among the 0.1% of the population
  S_random[1] <- sum(x[index2]==0);E_random[1] <- sum(x[index2]==1) 
  
  for (i in 2:nt) { # loop over days
    u <- runif(n) # uniform random deviates
    Iota<-sum(beta[x==2])#sum all currently I state people
    x[x==2&u<delta] <- 3 # I -> R with prob delta - 1/5
    x[x==1&u<gamma] <- 2 # E -> I with prob gamma - 1/3
    x[x==0&u<lambda*beta*Iota] <- 1 # S -> E with prob lambda*βj*sum(βi)
    
    # count the number of S state people at the ith day of 3 groups separately
    S[i] <- sum(x==0);S_low[i]<- sum(x[index1]==0);S_random[i]<- sum(x[index2]==0)
    
    # calculate the number of new infections at the ith day of 3 groups by computing the difference between the (i-1)th day and ith day
    E_new[i]<- S[i-1]-S[i];E_low[i]<- S_low[i-1]-S_low[i];E_random[i]<- S_random[i-1]-S_random[i]
  }
  list(E_new=E_new,E_low=E_low,E_random=E_random) 
} ## seir

set.seed(3)
ep <- seir() # run simulation 


# Step 2: Standardize
# standardize each trajectory by calculating the incidence per 10000 people per day
# Incidence is calculated by dividing their size
n=5.5e+6
ep_new <- ep$E_new/n*10000 ## whole population
ep_low <- ep$E_low/n/0.1*10000 ## 10% of the population
ep_random <- ep$E_random/n/0.001*10000 ## 0.1% of the population


# Step 3: Plot daily trajectory of new infections 
# plot and draw the trajectory of new infections with peak in a random sample of 0.1% of the population
plot(ep_random,ylim=c(0,max(ep_new,ep_low,ep_random)),# set up the y-axis limit of the graph
     main="Daily Infection Trajectories",xlab="day",ylab="Incidence per 10000 per day", # write the title of the plot and label axes
     type="l",col="dodgerblue",lwd=3) # using blue line to display the trajectory and assign the line width

# add guideline to fix the position of peak value; average the x-axis value in case of multiple peak values
abline(h=max(ep_random), v=mean(which(ep_random == max(ep_random))), lty=2, col="dodgerblue",lwd=2)

# mark the point of peak value of the trajectory
points(mean(which(ep_random == max(ep_random))),max(ep_random), col="brown",cex=1.5,pch=16)

# label the day on which peak occurs at the near place
text(mean(which(ep_random == max(ep_random)))-2,max(ep_random)-6, mean(which(ep_random == max(ep_random))),col="brown",cex=1.2) 

# Draw the trajectory of new infections with peak among the 10% of the population similarly
lines(ep_low,col="grey",lwd=4)
abline(h=max(ep_low), v=mean(which(ep_low == max(ep_low))), lty=2, col= "grey",lwd=2)
points(mean(which(ep_low == max(ep_low))),max(ep_low), col="brown",cex=1.5,pch=16)
text(mean(which(ep_low == max(ep_low)))+2,max(ep_low)+6, mean(which(ep_low == max(ep_low))),col="brown",cex=1.2)

# Draw the trajectory of new infections with peak among the whole population as above
lines(ep_new, col="black",lwd=3)
abline(h=max(ep_new), v=mean(which(ep_new == max(ep_new))), lty=2, col="black", lwd=2)
points(mean(which(ep_new == max(ep_new))),max(ep_new), col="brown",cex=1.5,pch=16)
text(mean(which(ep_new == max(ep_new)))+2,max(ep_new)+6, mean(which(ep_new == max(ep_new))),col="brown",cex=1.2)
 

# classify different groups by adding legend
legend(list(x=2,y=170),legend=c("0.1% random","cautious 10%","whole population","day of peak"),
       col=c("dodgerblue","grey","black","brown"),lty=1,pch=c(NA,NA,NA,16),lwd=2,cex=0.9,bty="n")


# Step 4: 10 repeated simulations
# store the corresponding infection and peak values at each day in list
new <- low <- random <- list() 
peak_whole <- peak_low <- peak_random <-rep(0,10)
for (i in 1:10) {
  n=5.5e+6
  ep<-seir()
  new[[i]] <- ep$E_new/n*10000 
  low[[i]] <- ep$E_low/n/0.1*10000
  random[[i]]<- ep$E_random/n/0.001*10000
  peak_whole[i]<-max(new[[i]])
  peak_low[i]<-max(low[[i]])
  peak_random[i]<-max(random[[i]])
}


# Step 5: Visualize(zoom)
# There are 3 plots showing the results of 10 repeated simulations of each group respectively.
# plot and draw the first trajectory for new infections each day among the whole population
par(mfrow=c(2,2), mar=c(4,4,1,1),oma = c(0, 0, 2, 0)) # set margin and outer margin of the graph and create 4 regional graphs inside
plot(new[[1]],xlab="day",ylab="Incidence per 10000 per day", 
     type="l",col="black",lwd=2) 
legend("topleft","Whole Population",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")# add legend
max_new = mean(which(new[[1]]==max(new[[1]])))# Find and store the day of first peak value occurs
for (i in 2:10) {# loop
  lines(new[[i]],col="black",lwd=2)# draw the other 9 new infections among whole population
  max_new  = max_new + mean(which(new[[i]]==max(new[[i]]))) #sum the day on which peak values occurs of 10 simulations
}
abline(v=max_new/10, lty=2, col="red", lwd=2)# average the day of total value and add a guideline to indicate
axis(side=3, at=max_new/10,labels=max_new/10,col ="red", col.axis="red")# label the average value by axis

# draw trajectory for a random sample of 0.1% of the population as above
plot(random[[1]],xlab="day",ylab="Incidence per 10000 per day", type="l",col="dodgerblue",lwd=2) 
legend("topleft","0.1% Random",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")
max_random = mean(which(random[[1]]==max(random[[1]])))
for (i in 2:10) {
  lines(random[[i]],col="dodgerblue",lwd=2)
  max_random = max_random + mean(which(random[[i]]==max(random[[i]])))
}
abline(v=max_random/10, lty=2, col="red", lwd=2)
axis(side=3, at=max_random/10,labels=max_random/10,col="red",col.axis="red")

# draw trajectory for the 10% of the population similarly
plot(low[[1]],xlab="day",ylab="Incidence per 10000 per day", type="l",col="grey",lwd=2,) 
legend("topleft","Cautious 10%",text.col="red",lty=0,lwd=2,cex=0.9,bty ="n")
max_low = mean(which(low[[1]]==max(low[[1]])))
for (i in 2:10) {
  lines(low[[i]],col="grey",lwd=2)
  max_low = max_low + mean(which(low[[i]]==max(low[[i]])))
}
abline(v=max_low/10, lty=2, col="red", lwd=2)
axis(side=3, at=max_low/10,labels=max_low/10,col="red",col.axis="red")

# There are 3 plots showing trajectories of peak values of 10 repeated simulations
# Plot and draw the first trajectory for peak values of new infections each day among the whole population
plot(peak_whole,ylim=c(0,max(peak_whole,peak_low,peak_random)), # set up the y-axis limit of the graph and label axes
     xlab="simulation",ylab="peak incidence per 10000 per day", type="l",col="black",lwd=2)
# Using 'lines' to draw trajectory for the 10% of the population & a random sample of 0.1% of the population
lines(peak_low,col="grey",lwd=2); lines(peak_random,col="dodgerblue",lwd=2);
legend("bottom",legend=c("0.1% random","cautious 10%","whole population"),
       bty="n",col=c("dodgerblue","grey","black"),lty=1,lwd=2,cex=0.5, horiz=TRUE)# add legend 
mtext("Variability of 10 Repeated Simulations", side = 3, line = 0, outer = T)# add main title of the whole plot


# Step 6: Comment
# People who choose to download ZOE app are more cautious about Covid than average. The ZOE sample is with lower transmission probabilities. 
# The general trend of the Incidence of ZOE sample and whole population are similar: the incidence rates increase at first until they reach the peak value and then decrease.
# However, peak values and the days on which peak values occur are different: the ZOE sample's peak value day will occur later with generally smaller peak value than the whole population.
