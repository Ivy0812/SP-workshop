

bfgs <- function (theta,f,...,tol=1e-5,fscale=1,maxit=100){
  
  first <- function(theta,f,...){
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
  
  if (is.null(attributes(f(theta,...)))){ 
    fd <- first(theta,f,...)
  }   
  else{
    fd <- attr(f(theta,...),"gradient")
  }
      
  iter <- 0
  Bk1 <- diag(length(theta))
  for (i in 1:maxit){
    sk <- delta <- - Bk1 %*% first(theta,f,...)
    theta1 <- theta + delta
    yk <- first(theta1,...) - first(theta,f,...)
    pk <- 1/t(sk) %*% yk
    Bk1 <- (diag(length(theta))- pk*sk %*% t(yk)) %*% Bk1 %*% (diag(len(th0))- pk*yk %*% t(sk)) + pk*sk %*% t(sk)
    iter <- iter+1
    if (max(abs(first(theta,f,...))) < (abs(f(theta,...))+fscale)*tol) {
      stop()
    }
    else{
      theta <- theta1
    }
  } 
  gll0 <- first(theta,f,...) ## gran of nll at th0 
  eps <- 1e-7 ## finite difference interval
  th0 <- theta
  Hfd <- matrix(0,length(th0),length(th0)) ## finite diference Hessian 
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps 
    gll1 <- first(th1,f,...) ## compute resulting nll
    Hfd[i,] <- (gll1 - gll0)/eps ## approximate second derivs
      }     
  
   list(f=f(theta,...), theta=theta, iter=iter, g=first(theta,f,...), H=Hfd)     
  }
    


