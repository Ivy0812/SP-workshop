

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  
  if (!(is.finite(f(theta,...))))  stop("   ")
  
  
  if (is.null(attributes(f(theta,...), "gradient"))){
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
  } 
  else{
    first <- function(theta,f,...) {
      fd <- attr(f(theta,...),"gradient")
      fd
    }
  }
  
  
  # 1. If the objective or derivatives are not finite at the initial theta
  for (i in first(theta,f,...)){
    if (!(is.finite(i))) stop("   ")
  }
 
  
  Bk1 <- diag(length(theta))
  for (iter in 1:maxit){
    if (max(abs(first(theta,f,...))) < (abs(f(theta,...))+fscale)*tol) break
    else{
      delta <- - Bk1 %*% first(theta,f,...)
      while (!(is.finite(f(theta+delta,...))) | (f(theta,...) < f(theta+delta,...))) delta=delta/2
      # 2. If the step fails to reduce the objective but convergence has not occurred
      k <- 0.1
      while (t(first(theta+delta,f,...)) %*% delta < 0.9 * t(first(theta,f,...)) %*% delta){
        delta = 1 + k * delta
        if (!(is.finite(f(theta+delta,...))) | (f(theta,...) < f(theta1,...))) {
            k <- k / 2
            if (k < 0.001) stop("   ")
        }
      }
        
      yk <- first(theta + delta,f,...) - first(theta,f,...)
      pk <- drop(1/t(delta) %*% yk)
      A <- Bk1 - pk*delta %*% (t(yk) %*% Bk1)
      Bk1 <- A - pk*(A %*% yk) %*% t(delta) + pk*delta %*% t(delta)
      theta <- theta + delta
    }
  }

  # 3. If maxit is reached without convergence. �?
  if ((max(abs(first(theta,f,...))) >= (abs(f(theta,...))+fscale)*tol)) stop(" ")
  
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









#20211118
# This work is completed by group 12：Wenyi Yu s2161093, Minke Pan s2160782 and Yanren Mao s2207399


# Create a function bfgs with input theta, f, ..., tol=1e-5, fscale=1, maxit=100
# f is an objective function with arguments: 
# the first argument is the vector of optimization parameters;
# second argument of f is a logical of whether the gradient of f should be computed
# remaining arguments are passed by "..."
# The function bfgs will implement a quasi-Newton minimization method and return a named list containing
# 1. f the scalar value of the objective function at the minimum
# 2. theta the vector of values of the parameters at the minimum
# 3. iter the number of iterations taken to reach the minimum.
# 4. g the gradient vector at the minimum (so the user can judge closeness to numerical zero).
# 5. H the approximate Hessian matrix (obtained by finite differencing) at the minimum

###################################################################################################



# Step 1: 
# Firstly, check if the objective function f is infinite:
##if it is infinite, then stop 
# Secondly, check if f has a gradient atrribute:
##If no, then we are going to implement the gradient ourselves using finite approximations; 
##If yes, simply use "attr" can help us to obtain the gradient

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){
  
  if (!(is.finite(f(theta,...))))  stop("the objective function is not finite")  # check if the objective function is infinite, if true then stop and issue error
  
  
  if (is.null(attr(f(theta,...), "gradient"))){ # check if f has null attribute, if true then continute to compute gradient
    # Finite difference approximations
    first <- function(theta,f,...) {
      fd <- th0 <- theta # test parameter value and initial theta value
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
  } 
  else{# if the function has gradient attribute
    first <- function(theta,f,...) {
      fd <- attr(f(theta,...),"gradient") # use "attr" to obtain the gradient attribute
      fd
    }
  }

  ###################################################################################################



  # Step 2:
  # Now we need to issue errors for bfgs function at the following cases:
  # 1. If the objective or derivatives are not finite at the initial theta
  ##use "is.finite" to check
  # 2. If the step fails to reduce the objective but convergence has not occurred
  # 3. If maxit is reached without convergence

  # and also need to consider:
  # 1. need to reduce the step length if a step leads to a non-finite objective value and 
  ##also need to make sure step length reduces the objective function and satisfies the second wolfe condition(c_2 = 0.9)
  # 2. consider the magnitude of the objective
  ##suitable convergence -- max(abs(g)) < (abs(f0)+fscale)*tol where g is the current gradient and f0 is the current value of the objective
  # 3. the final approximate Hessian matrix can be generated by finite differencing the gradient vector
  ##need to be asymmetric -- H <- 0.5 * (t(H) + H) 
  # 4. the cost of the approximate inverse Hessian should be O(p^2)
  # 怎么来的？？？
  # we will include the above error issuing and considerations in the following 

###################################################################################################



  # Step 2.1
  # In the following section, we are going to 
  # A. issue the error when the objective or derivatives are not finite at the initial theta
  # B. consideration 1:
    ## the magnitude of f to obtain a suitable convergence
    ## reduce the step length if a step leads to a non-finite objective value
    ## step length reduces the objective function and satisfies the second wolfe condition(c_2 = 0.9)
  ### a. If the step length fails to reduce objective function or the objective is non-finite, need to reduce step length
  ### b. If reducing the step length in a. is unlikely to meet Wolfe condition, this means step length need to be increased slightly *(1+k) where k is the increment
  ### c. If increasing the step length in b. fails to reduce the objective or objective is non-finite, need to halve the step length
  ### d. If no amount of reduction gives a reduction of finite objective, need to stop and issue an error
  
 
  # If the objective or derivatives are not finite at the initial theta
  for (i in first(theta,f,...)){ 
    if (!(is.finite(i))) stop("the objective or derivatives are not finite at the initial theta")
  }
 
  # consider the magnitude of the objective
  # suitable condition for convergence is max(abs(g)) < (abs(f0)+fscale)*tol 
  # subsititute theta in 'first' function to get the current gradient value
  # subsititute theta in f to get the current objective value 
  Bk1 <- diag(length(theta))
  for (iter in 1:maxit){ # loop over max iterations
    if (max(abs(first(theta,f,...))) < (abs(f(theta,...))+fscale)*tol) break #  step length both reduces the objective function with suitable convergenece then break
    else{# if convergence does not occur
      delta <- - Bk1 %*% first(theta,f,...) # the initial Quasi-Newton step is ∆ = −B[k]∇D(θ[k])
      # 2. If the step fails to reduce the objective but convergence has not occurred
      while (!(is.finite(f(theta+delta,...))) | (f(theta,...) < f(theta+delta,...))) delta <- delta/2 # if the objective is non-finite or fails to reduce the objective

      k <- 0.1 #increment of step length
      # If the second Wolfe condition is not met, i.e. ∇D(θ[k] + ∆)^T∆ < c2∇D(θ[k])^T∆ 
       ##(the second Wolfe condition ∇D(θ[k] + ∆)^T∆ ≥ c2∇D(θ[k])^T∆ ensures B[K+1] stays positive definite)
      while (t(first(theta+delta,f,...)) %*% delta < 0.9 * t(first(theta,f,...)) %*% delta){
        delta <- (1 + k) * delta # increase step length by k
        # if reducing the step length leads to non-finite obejective or fails to reduce the objective, need to reduce step length
        if (!(is.finite(f(theta+delta,...))) | (f(theta,...) < f(theta+delta,...))) {
            k <- k / 2 # so need to decrease the step length by halving it
            if (k < 0.001) stop("the step fails to reduce the objective") # If no amount of reduction gives a reduction of finite objective, stop and issue an error
        }
      }

      #Step 2.2 -- cost问minke
      yk <- first(theta + delta,f,...) - first(theta,f,...)
      pk <- drop(1/t(delta) %*% yk)
      A <- Bk1 - pk*delta %*% (t(yk) %*% Bk1)
      Bk1 <- A - pk*(A %*% yk) %*% t(delta) + pk*delta %*% t(delta)
      theta <- theta + delta
    }
  }



  # Step 2.3:
  # In the following section, we are going to issue error when:
  # maxit is reached without convergence - will jump out of the loop
  # There are two cases of the second for loop(!label line # in final version) in section 2.1
  # a. if convergence occurs , then break：jump out of the for loop and will not run the following code to issue error
  # b. if convergence does not occur, so continue the process in for loop, when reach maxit but still no convergence，jump out of for loop and run the following code
  if ((max(abs(first(theta,f,...))) >= (abs(f(theta,...))+fscale)*tol)) stop("maxit is reached without convergence")
  


  # Step 2.4:
  # In the following section, we mainly do 2 things:
  # A. approximate Hessian matrix (obtained by finite differencing) at the minimum
  # B. make sure Hessian matrix meet inequality H <- 0.5 * (t(H) + H) to ensure asymmetricity
  gll0 <- first(theta,f,...) ## gran of nll at th0 
  eps <- 1e-7 ## finite difference interval
  th0 <- theta
  Hfd <- matrix(0,length(th0),length(th0)) ## finite difference Hessian 
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
