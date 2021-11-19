# This work is completed by group 12：Wenyi Yu s2161093, Minke Pan s2160782 and Yanren Mao s2207399
###################################################################################################

# Overview

# The following code:

# 1.create a function "bfgs" with the input of: theta, f, ..., tol=1e-5, fscale=1, maxit=100
# f is an objective function with arguments: 
# a) the first argument is the vector of optimization parameters
# b) second argument of f is a logical of whether the gradient of f should be computed
# c) remaining arguments are passed by "..."

# 2.if the supplied f does not supply a gradient attribute when requested, then compute the gradient 
# by finite differencing. 

# 3.issue errors or warnings (using stop or warning as appropriate) in at the following cases:
# a) If the objective or derivatives are not finite at the initial theta
# b) If the step fails to reduce the objective but convergence has not occurred 
# c) If maxit is reached without convergence

# 4.implement a quasi-Newton minimization method and return a named list containing:
# a) f - the scalar value of the objective function at the minimum
# b) theta - the vector of values of the parameters at the minimum
# c) iter - the number of iterations taken to reach the minimum
# d) g - the gradient vector at the minimum (so the user can judge closeness to numerical zero)
# e) H - the approximate Hessian matrix (obtained by finite differencing) at the minimum
########################################################################################

# Step 1: compute gradient if needed

# Firstly, check if the objective function f is infinite; if so, stop the function
# Secondly, check if f has a gradient attribute:
##If no, create a function "first" using finite approximations to calculate the gradient; 
##If yes, simply use "attr" to obtain the gradient and put it into the independent function, 
# which is also named "first" taking the following loop into account.
####################################################################

bfgs <- function(theta,f,...,tol=1e-5,fscale=1,maxit=100){ # create a function named "bfgs"
  
  # if the objective value at the initial theta is infinite, "is.finite" will return FALSE
  # "!FALSE" means TRUE, then stop and issue error
  if (!(is.finite(f(theta,...))))  stop("the objective value is not finite") 
  
  # if the supplied f does not supply a gradient attribute when requested, 
  # "is.null" will return TURE then compute gradient
  if (is.null(attr(f(theta,...), "gradient"))){ 
    # create a function named "first" for finite difference approximations
    first <- function(theta,f,...) {
      fd <- th0 <- theta # test parameter value and initial theta value
      nll0 <- f(th0,...) ## nll at th0
      eps <- 1e-7 ## finite difference interval
      for (i in 1:length(th0)) { ## loop over parameters
        th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps 
        nll1 <- f(th1,...) ## compute resulting nll
        fd[i] <- (nll1 - nll0)/eps ## approximate -dl/dth[i]
      }
      fd
    }
  } 
  else{# if the supplied f has a gradient attribute
    # use "attr" to obtain the gradient value and put it into independent function "first"
    first <- function(theta,f,...) attr(f(theta,...),"gradient")
  }
################################################################
  
# Step 2: check some errors and adjust step length

# errors should be issued for bfgs function at the following cases:
# 1. if the derivatives are not finite at the initial theta
# 2. if the step leads to a non-finite objective value or fails to reduce the objective 
# but convergence has not occurred

# consider the magnitude of the objective as well as suitable convergence by the code:
# max(abs(g)) < (abs(f0)+fscale)*tol 
# where g is the current gradient and f0 is the current value of the objective 
  
# step length should be adjusted:
# 1. reduce the step length if a step leads to a non-finite objective value or fails to reduce the objective
# 2. increase step length if it does not meet the second wolfe condition (c_2 = 0.9)
# 3. if the increased step length leads to a non-finite objective value or fails to reduce the objective again,
# increase the step length less and check again
# 4 repeat the above steps, if no amount of reduction gives a reduction, stop and issue a warning.
##################################################################################################
  
  # if the derivatives are not finite at the initial theta, "is.finite" will return FALSE
  # "!FALSE" means TRUE, then stop and issue error
  for (i in first(theta,f,...)){ 
    if (!(is.finite(i))) stop("at least 1 derivative is not finite at the initial theta")
  }
  
  # set the initial inverse approximate Hessian to Identity matrix to avoid the need for a matrix solve
  # with row and column numbers equal to the length of theta
  Bk <- diag(length(theta))
  for (iter in 1:maxit){ # loop over the maximum number of iterations before giving up
    # consider the magnitude of the objective and if the step length reduces the objective function with 
    # suitable convergence, break and jump to line 160
    if (max(abs(first(theta,f,...))) < (abs(f(theta,...))+fscale)*tol) break
    else{# if convergence does not occur
      delta <- - Bk %*% first(theta,f,...) # compute the initial Quasi-Newton step - ∆ = −B[k]∇D(θ[k])
      # if the step leads to a non-finite objective value or fails to reduce the objective 
      # repeatedly halve ∆ until D(θ + ∆) < D(θ)
      while (!(is.finite(f(theta+delta,...))) | (f(theta,...) < f(theta+delta,...))) delta <- delta/2 

      k <- 0.1 # set initial increment of step length to 0.1
      step <- delta # fix delta and change step by "step"
      # if the second Wolfe condition is not met, i.e. ∇D(θ[k] + ∆)^T∆ < c2∇D(θ[k])^T∆ 
      ##(the second Wolfe condition ∇D(θ[k] + ∆)^T∆ ≥ c2∇D(θ[k])^T∆ ensures B[K+1] stays positive definite)
      while ((t(first(theta+step,f,...)) %*% step) < (0.9 * t(first(theta,f,...)) %*% step)){
        step <- (1 + k) * delta # increase step length by k * delta
        # if reducing the step length leads to non-finite objective or fails to reduce the objective, 
        # reduce the increment k
        if (!(is.finite(f(theta+step,...))) | (f(theta,...) < f(theta+step,...))) {
            k <- k / 2 # repeatedly halve k 
            # if no amount of reduction gives a reduction of finite objective, the while will loop forever
            # so set the floor constraint, once k < 0.001, stop and issue an error
            if (k < 0.001) stop("fail to find out suitable step length for convergence after adjustments") 
        }
      }
###########################################################################################################
      
# Step 3: updates the inverse approximate Hessian - Quasi-Newton

# 1.Let θ[k] denote the kth trial θ, with approx. inverse Hessian B[k].
# 2.Let s[k] = θ[k+1] − θ[k] and y[k] = ∇D(θ[k+1]) − ∇D(θ[k])
# 3.Defining ρ[k]^(-1) = s[k]^T*y[k] the BFGS update is as following (since ρk−1 is a constant, ρk can be simply calculated by reciprocal)
# 4.The Quasi-Newton step from θ[k] is ∆ = −B[k]∇D(θ[k]), which may adjusted by Step 2  
# update of the approximate inverse Hessian to have cost O(p2) (p - the number of parameters) 
# 5.B[k+1] = (I-ρ[k]s[k]y[k]^T)B[k](I-ρ[k]y[k]s[k]^T)+ρ[k]s[k]s[k]^T, expand bracket: 
#          = (B[k]-ρ[k]s[k]y[k]^T*B[k])(I-ρ[k]y[k]s[k]^T) + ρ[k]s[k]s[k]^T
#   Let A  = B[k]-ρ[k]s[k]y[k]^T*B[k] then
#   B[k+1] = A(I-ρ[k]y[k]s[k]^T) + ρ[k]s[k]s[k]^T
#          = A - ρ[k](Ay[k])(s[k]^T) + ρ[k]s[k]s[k]^T 
#####################################################  
      
      delta <- step # assign step value to delta variable
      yk <- first(theta + delta,f,...) - first(theta,f,...) # y[k] = ∇D(θ[k+1]) − ∇D(θ[k])
      pk <- drop(1/t(delta) %*% yk) # since p[k] is a 1*1 matrix, need to use "drop" to extract constant value
      A <- Bk - pk*delta %*% (t(yk) %*% Bk) # define A = B[k]-ρ[k]s[k]y[k]^T*B[k]
      Bk <- A - pk*(A %*% yk) %*% t(delta) + pk*delta %*% t(delta) # expand the bracket and substitute A 
      theta <- theta + delta # update theta
    }
  }
#############################

#Step 4: generate final approximate Hessian matrix

# 1.before calculation, check if maxit is reached without convergence
# 2.generate final approximate Hessian matrix by finite differencing the gradient vector.  
# 3.ensure approximate Hessian will be asymmetric
# 4.output the required (named) list
####################################
  
  # if maxit is reached without convergence - issue error and jump out of the loop
  # There are two cases of the iterations in line 95 in Step 2
  # a. if convergence occurs, the condition of break is opposite to the following condition in line 160 
  # so function will not be stopped here.
  # b. if convergence does not occur, for loop will end by defaut twhen maxit is reached，which means 
  # the following condition is satisfied, so function will be stopped.
  if ((max(abs(first(theta,f,...))) >= (abs(f(theta,...))+fscale)*tol)) stop("maxit is reached without convergence")
  
  # approximate Hessian matrix (obtained by finite differencing) at the minimum
  gll0 <- first(theta,f,...) ## gran of nll at th0 
  eps <- 1e-7 ## finite difference interval
  th0 <- theta
  Hfd <- matrix(0,length(th0),length(th0)) ## finite difference Hessian 
  for (i in 1:length(th0)) { ## loop over parameters
    th1 <- th0; th1[i] <- th1[i] + eps ## increase th0[i] by eps 
    gll1 <- first(th1,f,...) ## compute resulting nll
    Hfd[i,] <- (gll1 - gll0)/eps ## approximate second derivs
  }     
  H <- 0.5 * (t(Hfd) + Hfd) # make sure Hessian matrix is symmetric
  
  # if the supplied f has a gradient attribute, hide the attribute when output
  D <- f(theta,...) # set D to the final minimum objective value
  attr(D,"gradient") <- NULL # set the "gradient" attribute of D to NULL
  # return the required list, iter need to reduce by 1 since the loop starts by iter = 1 without actual iteration
  list(f=D, theta=theta, iter=iter-1, g=first(theta,f,...), H=H)     
}
