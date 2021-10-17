seir <- function(n=5.5e+6,ne=10,nt=150,gamma=1/3,delta=1/5) {
  ## SEIR stochastic simulation model.
  ## n = population size; ni = initially infective; nt = number of days ## gamma = daily prob E -> I; delta = daily prob I -> R;
  ## bmu = mean beta; bsc = var(beta) = bmu * bsc
  lambda <- 0.4/n
  x <- rep(0,n) ## initialize to susceptible state
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta) ## individual infection rates 
  index1 <- order(beta)[1:(0.1*n)]
  index2 <- sample(n,0.01*n,replace=FALSE)
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
ep <- seir() ## run simulation 
ep_new <- (ep$E_new-min(ep$E_new))/(max(ep$E_new)-min(ep$E_new))*100
ep_low <- (ep$E_low-min(ep$E_low))/(max(ep$E_low)-min(ep$E_low))*100
ep_random <- (ep$E_random-min(ep$E_random))/(max(ep$E_random)-min(ep$E_random))*100

plot(ep_new,ylim=c(0,100),xlab="day",ylab="Incidence", type="l",lwd=2) ## E black
lines(ep_low,col=4,lwd=2); lines(ep_random, col=2,lwd=2)

  