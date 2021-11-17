

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  if (is.infinite(f(theta,...)))  stop("   ")
  if (is.null(attributes(f(theta,...)))){
    first <- function(theta,f,...) {
      fd <- th0 <- theta
      nll0 <- f(theta,...) ## nll at th0
      eps <- 1e-7 ## finite difference interval
      for (i in 1:length(th0)) { ## loop over parameters
        th1 <- th0; th1[i] <- th1[i] + eps
        th1[i] <- th1[i] + eps ## increase th0[i] by eps 
        nll1 <- f(th1,...) ## compute resulting nll
        fd[i] <- (nll1 - nll0)/eps ## approximate -dl/dth[i]
      }
    fd
    }
  } else{
    first <- function(theta,f,...) {
    fd <- attr(f(theta,...),"gradient")
    fd
    }
  }
  # 1. If the objective or derivatives are not finite at the initial theta
  for (i in first(theta,f,...)){
    if (is.infinite(i)) stop("   ")
  }
 
  Bk1 <- diag(length(theta))
  eps <- 1e-7
  for (iter in 1:maxit){
    sk <- delta <- - Bk1 %*% first(theta,f,...)
    theta1 <- theta + delta
    yk <- first(theta1,f,...) - first(theta,f,...)
    pk <- drop(1/t(sk) %*% yk)
    A <- Bk1 - pk*sk %*% (t(yk) %*% Bk1)
    Bk1 <- A - pk*(A %*% yk) %*% t(sk) + pk*sk %*% t(sk)

    if (max(abs(first(theta,f,...))) < (abs(f(theta,...))+fscale)*tol) {
      break
    } else{
      # 2. If the step fails to reduce the objective but convergence has not occurred
      if (f(theta,...) < f(theta1,...)) stop("")
      theta <- theta1
    }
}

  # 3. If maxit is reached without convergence. 丑
  if ((iter == maxit) & (max(abs(first(theta,f,...))) >= (abs(f(theta,...))+fscale)*tol)) stop(" ")
  
  gll0 <- first(theta,f,...) ## gran of nll at th0 
  eps <- 1e-7 ## finite difference interval
  th0 <- theta
  Hfd <- matrix(0,length(th0),length(th0)) ## finite diference Hessian 
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps 
    gll1 <- first(th1,f,...) ## compute resulting nll
    Hfd[i,] <- (gll1 - gll0)/eps ## approximate second derivs
  }     
  H <- 0.5 * (t(Hfd) + Hfd)
  list(f=f(theta,...), theta=theta, iter=iter, g=first(theta,f,...), H=Hfd)     
}

rb <- function(theta,getg=FALSE,k=10) {
  ## Rosenbrock objective function, suitable for use by ’bfgs???
  z <- theta[1]; x <- theta[2]
  f <- k*(z-x^2)^2 + (1-x)^2 + 1
  if (getg) {
    attr(f,"gradient") <- c(2*k*(z-x^2),
                            -4*k*x*(z-x^2) -2*(1-x))
  }
  f
} ## rb

theta<-c(-1,2)
bfgs(theta,rb,tol=1e-5,fscale=1,maxit=100) 

#debug(bfgs(theta,rb,getg=FALSE,k=10,tol=1e-5,fscale=1,maxit=100))
