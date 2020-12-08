
###############################################
## Functions Needed for Parameter Estimation ##
###############################################


## Load the library necessary for non-linear optimization
library(Rsolnp)

## The following function computes the maximum likelihood estimate of the 
## HDCSBM parameters. Note that it depends on the functions theta.con() 
## and log.likelihood() defined below.
get.MLE <- function(Adj, n.comm, membership){
  #################################################################################################
  #  This function computes the maximum likelihood estimates associated with the HDCSBM           #
  #                                                                                               #
  # Inputs:                                                                                       #
  #         Adj = the nxn (symmetric, undirected) adjacency matrix                                #
  #      n.comm = an array of length B whose elements define the size of each                     #
  #               of the B communities                                                            #
  #  membership = an array of length n whose elements are the numbers 1,2,...,B                   #
  #               indicating which of the B communities each node belongs to                      #
  #                                                                                               #
  # Outputs:                                                                                      #
  #   estimates = an array of length B(B+1)+n of pi estimates followed by lambda estimates        #
  #               followed by theta estimates. In the case of the pi's and lambda's , the         #
  #               results are ordered by block pairs as follows: (1,1), (1,2), ..., (1,B),        #
  #               (2,2), (2,3), ..., (2,B), ..., (3,3), (3,4), ..., (3,B), ..., (B,B)             #
  #################################################################################################
  # Get the data-related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  Delta <- (Adj == 0) * 1 # delta matrix indicating which elements of the adjacency matrix are 0
  
  # Estimates of the pis:
  pi.hat <- rep(0, num.blocks)
  block.count <- 1
  for(i in 1:B){
    for(j in i:B){
      indx.i <- which(membership == i)
      indx.j <- which(membership == j)
      if(i == j){
        pi.hat[block.count] <- mean(Delta[indx.i, indx.j][upper.tri(Delta[indx.i, indx.j])])
      }else{
        pi.hat[block.count] <- mean(Delta[indx.i, indx.j])
      }
      block.count <- block.count + 1
    }
  }
  all.zero <- which(pi.hat == 1)
  
  # Estimates of the lambdas and thetas:
  res <- solnp(pars = c(rep(4, num.blocks), rep(1, n)), 
               fun = log.likelihood, eqfun = theta.con, 
               eqB = rep(0, B), LB = rep(0.01, num.blocks+n), UB = rep(Inf, num.blocks+n), 
               control=list(trace=0), 
               Adj = Adj, n.comm = n.comm, membership = membership)
  lambda.hat <- res$par[1:num.blocks]
  lambda.hat[all.zero] <- NA
  theta.hat <- res$par[(num.blocks+1):(num.blocks+n)]
  
  # Return all estimates to the user
  return(list(estimates = c(pi.hat, lambda.hat, theta.hat)))
}



## The following functionalizes the constraint on the thetas
## (the thetas in a given community sum to that communitiy's size)
theta.con <- function(param, Adj, n.comm, membership){
  #################################################################################################
  # Inputs:                                                                                       #
  #      param = a parameter array of length B(B+1)/2 + n, where the first B(B+1)/2 entries       #
  #              correspond to the lambda's and the final n entries correspond to the theta's).   #
  #              Note that the lambda's should be ordered in block pairs as follows: (1,1),       #
  #              (1,2), ..., (1,B), (2,2), (2,3), ..., (2,B), ..., (3,3), (3,4), ..., (3,B),      #
  #              ..., (B,B)                                                                       #
  #         Adj = the nxn (symmetric, undirected) adjacency matrix                                #
  #      n.comm = an array of length B whose elements define the size of each                     #
  #               of the B communities                                                            #
  #  membership = an array of length n whose elements are the numbers 1,2,...,B                   #
  #               indicating which of the B communities each node belongs to                      #
  #                                                                                               #
  # Outputs:                                                                                      #
  #  constraint = the Lagrangian constraint function for the thetas                               #
  #################################################################################################
  # Get the data-related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  # Get the relevant parameters
  theta <- param[(num.blocks+1):length(param)]
  # Set up the constraint
  constraint <- rep(0,B)
  for(i in 1:B){
    constraint[i] <- sum(theta[which(membership == i)])-n.comm[i]
  }
  return(constraint) 
} 



## The following function defines the log-likelihood function 
## for the HDCSBM
log.likelihood <- function(param,  Adj, n.comm, membership){
  #################################################################################################
  #  This function computes the log-likelihood function associated with the hurdle DCSBM.         #
  #  Note that this ignores the pi's since they can be estimated separately.                      #
  #                                                                                               #
  # Inputs:                                                                                       #
  #      param = a parameter array of length B(B+1)/2 + n, where the first B(B+1)/2 entries       #
  #              correspond to the lambda's and the final n entries correspond to the theta's).   #
  #              Note that the lambda's should be ordered in block pairs as follows: (1,1),       #
  #              (1,2), ..., (1,B), (2,2), (2,3), ..., (2,B), ..., (3,3), (3,4), ..., (3,B),      #
  #              ..., (B,B)                                                                       #
  #         Adj = the nxn (symmetric, undirected) adjacency matrix                                #
  #      n.comm = an array of length B whose elements define the size of each                     #
  #               of the B communities                                                            #
  #  membership = an array of length n whose elements are the numbers 1,2,...,B                   #
  #               indicating which of the B communities each node belongs to                      #
  #                                                                                               #
  # Outputs:                                                                                      #
  #         -l = the negative of the log-likelihood function (we return a negative because the    #
  #              negative of the log-likelihood function is minimized)                            #
  #################################################################################################
  # Get the data-related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  A <- Adj
  Delta <- (A == 0) * 1 # delta matrix indicating which elements of the adjacency matrix are 0
  
  # Get the parameters 
  # BxB matrix for lambda (used to create the nxn one)
  lambda.mat <- matrix(0, nrow = B, ncol = B)
  lambda.mat[lower.tri(lambda.mat, diag = TRUE)] <- param[1:num.blocks]
  lambda.mat <- lambda.mat + t(lambda.mat)
  diag(lambda.mat) <- diag(lambda.mat)/2
  
  # nxn matrix for lambda
  Lambda <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      Lambda[i,j] <- lambda.mat[indx.r,indx.c]
    }
  }
  
  # Vector for thetas
  theta.vec <- matrix(param[(num.blocks+1):length(param)], nrow = n, ncol = 1)
  # nxn matrix of pairwise products of thetas
  Theta <- theta.vec %*% t(theta.vec)
  
  #  Matrix of ones
  I <- matrix(1, n, n)
  
  # Compute the log-likelihood:
  l <- (I-Delta)*(A*log(Lambda*Theta)-log(exp(Lambda*Theta)-I))
  l[lower.tri(l, diag = TRUE)] <- 0 # we need only the upper triangle
  l <- sum(l) # sum across the upper triangle
  
  return(-l) # return negative because we minimize the negative of the log-likelihood function
}
