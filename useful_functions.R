######################
## Useful Functions ##
######################

library(VGAM)
library(countreg)
library(Rsolnp)


# <BEGIN HDCSBM Adjacency Generator>
gen.Adj <- function(n.comm, membership, pi, lambda, theta){
  # number of nodes
  n <- sum(n.comm)
  #number of communities
  B <- length(n.comm)
  # little pi matrix
  pi.mat <- matrix(0, nrow = B, ncol = B)
  pi.mat[lower.tri(pi.mat, diag = TRUE)] <- pi
  pi.mat <- pi.mat + t(pi.mat)
  diag(pi.mat) <- diag(pi.mat)/2
  # big pi matrix 
  pi.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      pi.mat.BIG[i,j] <- pi.mat[indx.r,indx.c]
    }
  }
  # little lambda matrix
  lambda.mat <- matrix(0, nrow = B, ncol = B)
  lambda.mat[lower.tri(lambda.mat, diag = TRUE)] <- lambda
  lambda.mat <- lambda.mat + t(lambda.mat)
  diag(lambda.mat) <- diag(lambda.mat)/2
  # big lambda matrix 
  lambda.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      lambda.mat.BIG[i,j] <- lambda.mat[indx.r,indx.c]
    }
  }
  # matrix for theta * theta
  theta2.mat <- matrix(theta, ncol = 1) %*% matrix(theta, nrow = 1)
  
  # Adjacency matrix
  Adj <- rhpois(n = n^2, lambda = lambda.mat.BIG * theta2.mat, pi = 1-pi.mat.BIG)
  Adj[lower.tri(Adj)] <- 0
  Adj <- Adj + t(Adj)
  diag(Adj) <- 0
  
  return(Adj)
}
# </END HDCSBM Adjacency Generator>

# <BEGIN DCSBM Adjacency Generator>
gen.Adj.DCSBM <- function(n.comm, membership, lambda, theta){
  # number of nodes
  n <- sum(n.comm)
  #number of communities
  B <- length(n.comm)
  # little lambda matrix
  lambda.mat <- matrix(0, nrow = B, ncol = B)
  lambda.mat[lower.tri(lambda.mat, diag = TRUE)] <- lambda
  lambda.mat <- lambda.mat + t(lambda.mat)
  diag(lambda.mat) <- diag(lambda.mat)/2
  # big lambda matrix 
  lambda.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      lambda.mat.BIG[i,j] <- lambda.mat[indx.r,indx.c]
    }
  }
  # matrix for theta * theta
  theta2.mat <- matrix(theta, ncol = 1) %*% matrix(theta, nrow = 1)
  
  # Adjacency matrix
  Adj <- matrix(data = rpois(n^2, lambda.mat.BIG * theta2.mat), nrow = n, ncol = n)
  Adj[lower.tri(Adj)] <- 0
  Adj <- Adj + t(Adj)
  diag(Adj) <- 0
  return(Adj)
}
# </END DCSBM Adjacency Generator>


# <BEGIN Expected Adjacency>
exp.Adj <- function(n.comm, membership, pi, lambda, theta){
  # number of nodes
  n <- sum(n.comm)
  #number of communities
  B <- length(n.comm)
  # little pi matrix
  pi.mat <- matrix(0, nrow = B, ncol = B)
  pi.mat[lower.tri(pi.mat, diag = TRUE)] <- pi
  pi.mat <- pi.mat + t(pi.mat)
  diag(pi.mat) <- diag(pi.mat)/2
  # big pi matrix 
  pi.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      pi.mat.BIG[i,j] <- pi.mat[indx.r,indx.c]
    }
  }
  # little lambda matrix
  lambda.mat <- matrix(0, nrow = B, ncol = B)
  lambda.mat[lower.tri(lambda.mat, diag = TRUE)] <- lambda
  lambda.mat <- lambda.mat + t(lambda.mat)
  diag(lambda.mat) <- diag(lambda.mat)/2
  # big lambda matrix 
  lambda.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      lambda.mat.BIG[i,j] <- lambda.mat[indx.r,indx.c]
    }
  }
  # matrix for theta * theta
  theta2.mat <- matrix(theta, ncol = 1) %*% matrix(theta, nrow = 1)
  # nxn matrix of ones
  I <- matrix(data = 1, nrow = n, ncol = n)
  
  EA <- (I - pi.mat.BIG) * ((lambda.mat.BIG * theta2.mat) / (I - exp(-lambda.mat.BIG * theta2.mat)))
  diag(EA) <- 0
  EA[is.na(EA)] <- 0
  return(EA)
}
# </END Expected Adjacency>


# <BEGIN Log-Likelihood Function>
loglikelihood <- function(param,  Adj, n.comm, membership){
  #################################################################################################
  #  This function computes the log-likelihood function associated with the hurdle DCSBM.         #
  #  Note that this ignores the pi's since they can be estimated separately.                      #
  #                                                                                               #
  # Inputs:                                                                                       #
  #      param = a (B(B+1)/2 + n)x1 parameter vector (the first B(B+1)/2 entries are devoted to   #
  #              the lambda's and the final n entries are devoted to the theta's)                 #
  #        Adj = the nxn adjacency matrix                                                         #
  #     n.comm = a Bx1 vector whose elements indicate the size of each community                  #
  # membership = an nx1 vector indicating the community membership for each node                  #
  #                                                                                               #
  # Outputs:                                                                                      #
  #         -l = the negative of the gradients associated with each parameter (we return a        #
  #              negative because the negative of the log-likelihood function is minimized)       #
  #              This is a (B(B+1)/2 + n)x1 vector                                                #
  #################################################################################################
  # get data related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  A <- Adj
  Delta <- (A == 0) * 1 # delta matrix indicating which elements of the adjacency matrix are 0
  
  # get the parameters 
  # little matrix for lambda (used to create nxn ones)
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
  
  # vector for thetas
  theta.vec <- matrix(param[(num.blocks+1):length(param)], nrow = n, ncol = 1)
  # nxn matric of pairwise products of thetas
  Theta <- theta.vec %*% t(theta.vec)
  
  #  matrix of ones
  I <- matrix(1, n, n)
  
  # compute the log-likelihood:
  l <- (I-Delta)*(A*log(Lambda*Theta)-log(exp(Lambda*Theta)-I))
  l[lower.tri(l, diag = TRUE)] <- 0 # we need only the upper triangle
  l <- sum(l) # sum across the upper triangle
  
  return(-l) # return negative because we minimize the negative of the log-likelihood function
}
# </END Log-Likelihood Function>


# <BEGIN Theta Constraint>
theta.con <- function(param, Adj, n.comm, membership){
  # get data related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  # get relevant parameters
  theta <- param[(num.blocks+1):length(param)]
  constraint <- rep(0,B)
  for(i in 1:B){
    constraint[i] <- sum(theta[which(membership == i)])-n.comm[i]
  }
  return(constraint) 
} 
# </END Theta Constraint>


# <BEGIN ML Estimator>
get.MLE <- function(Adj, n.comm, membership){
  # get data related quantities
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
  
  # Standard errors of the pis:
  N <- n.comm %*% t(n.comm)
  diag(N) <- (diag(N) - n.comm)/2
  pi.mat <- matrix(0, nrow = B, ncol = B)
  pi.mat[lower.tri(pi.mat, diag = TRUE)] <- pi.hat
  pi.mat <- pi.mat + t(pi.mat)
  diag(pi.mat) <- diag(pi.mat)/2
  se.pi <- sqrt((pi.mat * (1-pi.mat)) / N)
  se.pi <- se.pi[lower.tri(se.pi, diag = TRUE)]

  # Estimates of the lambdas and thetas:
  initial.values <- c(rep(4, num.blocks), rep(1, n))
  res <- solnp(pars = initial.values, fun = loglikelihood, eqfun = theta.con, eqB = rep(0, B), LB = rep(0.01, num.blocks+n), UB = rep(Inf, num.blocks+n), control=list(trace=0), Adj = Adj, n.comm = n.comm, membership = membership)
  lambda.hat <- res$par[1:num.blocks]
  lambda.hat[all.zero] <- NA
  theta.hat <- res$par[(num.blocks+1):(num.blocks+n)]
  
  # # Standard errors of the lambdas and thetas:
  # hess <- res$hessian
  # se.lambda <- sqrt(diag(solve(hess)))[1:num.blocks]
  # se.lambda[all.zero] <- NA
  # se.theta <- sqrt(diag(solve(hess)))[(num.blocks+1):(num.blocks+n)]
  return(list(estimates = c(pi.hat, lambda.hat, theta.hat)))
  # return(list(estimates = c(pi.hat, lambda.hat, theta.hat), se = c(se.pi, se.lambda, se.theta)))
}
# </END ML Estimator>


# <BEGIN Gradient Function>
gradient <- function(param, Adj, n.comm, membership){
  #################################################################################################
  #   This function computes the gradients of the log-likelihood function associated with the     #
  #   hurdle DCSBM. Note that this ignores the pi's since they can be estimated separately.       #
  #                                                                                               #
  # Inputs:                                                                                       #
  #      param = a (B(B+1)/2 + n)x1 parameter vector (the first B(B+1)/2 entries are devoted to   #
  #              the lambda's and the final n entries are devoted to the theta's)                 #
  #        Adj = the nxn adjacency matrix                                                         #
  #     n.comm = a Bx1 vector whose elements indicate the size of each community                  #
  # membership = an nx1 vector indicating the community membership for each node                  #
  #                                                                                               #
  # Outputs:                                                                                      #
  #         -g = the negative of the gradients associated with each parameter (we return a        #
  #              negative because the negative of the log-likelihood function is minimized)       #
  #              This is a (B(B+1)/2 + n)x1 vector                                                #
  #################################################################################################
  # get data related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  delta.mat <- (Adj == 0) * 1 # delta matrix indicating which elements of the adjacency matrix are 0
  
  # get the parameters 
  # little matrix for lambda (used for blockwise iterations)
  lambda.mat <- matrix(0, nrow = B, ncol = B)
  lambda.mat[lower.tri(lambda.mat, diag = TRUE)] <- param[1:num.blocks]
  lambda.mat <- lambda.mat + t(lambda.mat)
  diag(lambda.mat) <- diag(lambda.mat)/2
  
  # big matrix for lambda (used for cellwise iterations)
  lambda.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      lambda.mat.BIG[i,j] <- lambda.mat[indx.r,indx.c]
    }
  }
  
  # vector for thetas (nx1)
  theta.vec <- param[(num.blocks+1):length(param)]
  
  # matrix for thetas (nxn theta1theta matrix)
  theta.mat <- matrix(theta.vec, ncol = 1) %*% matrix(theta.vec, nrow = 1)
  
  # pre-allocate vectors to store all gradients -- one for each parameter type
  grad.lambda.vec <- rep(0, num.blocks)
  grad.theta.vec <- rep(0, n)
  
  # gradients for lambda's
  block.count <- 1
  for(i in 1:B){
    for(j in i:B){
      indx.i = which(membership == i)
      indx.j = which(membership == j)
      # get block-specific values of lambda
      lambda <- lambda.mat[i,j]
      # gradients for lambda
      a <- Adj[indx.i, indx.j]
      delta <- delta.mat[indx.i, indx.j]
      theta1theta2 <- theta.mat[indx.i, indx.j]
      deriv.lambda <- (1-delta) * ((a/lambda) - theta1theta2/(1 - exp(-lambda*theta1theta2)))
      if(i == j){ # diagonal block gradients
        grad.lambda.vec[block.count] <- sum(deriv.lambda[upper.tri(deriv.lambda)])  
      }else{ # off-diagonal block gradients
        grad.lambda.vec[block.count] <- sum(deriv.lambda)  
      }
      block.count <- block.count + 1
    }
  }
  # gradients for theta's
  for(k in 1:n){
    theta_k <- theta.vec[k]
    # Adjacency T
    a <- matrix(0, nrow = n, ncol = n)
    a[k,] <- Adj[k,]
    a[,k] <- Adj[,k]
    # Delta T
    delta <- matrix(0, nrow = n, ncol = n)
    delta[k,] <- delta.mat[k,]
    delta[,k] <- delta.mat[,k]
    # Lambda T
    lambda <- matrix(0, nrow = n, ncol = n)
    lambda[k,] <- lambda.mat.BIG[k,]
    lambda[,k] <- lambda.mat.BIG[,k]
    # theta1theta2 T
    theta1theta2 <- matrix(0, nrow = n, ncol = n)
    theta1theta2[k,] <- theta.mat[k,]
    theta1theta2[,k] <- theta.mat[,k]
    # derivatives with respect  to theta_k
    deriv.theta <- (1-delta) * ((a/theta_k) - (lambda*(theta1theta2/theta_k))/(1-exp(-lambda*theta1theta2)))
    deriv.theta[is.nan(deriv.theta)] <- 0
    grad.theta.vec[k] <- 0.5*sum(deriv.theta) #divide by 2 here because we only need to sum the upper triangle
  }
  g <- c(grad.lambda.vec, grad.theta.vec) # full gradient vector
  return(-g) # return negative because we minimize the negative of the log-likelihood function
}
# </END Gradient Function>


# <BEGIN Observed Information Function>
obs.infomat <- function(param, Adj, n.comm, membership){
  #################################################################################################
  #   This function computes the observed information matrix with the hurdle DCSBM. Note that     #
  #   this ignores the pi's since they can be estimated separately.                               #
  #                                                                                               #
  # Inputs:                                                                                       #
  #      param = a (B(B+1)/2 + n)x1 parameter vector (the first B(B+1)/2 entries are devoted to   #
  #              the lambda's and the final n entries are devoted to the theta's)                 #
  #        Adj = the nxn adjacency matrix                                                         #
  #     n.comm = a Bx1 vector whose elements indicate the size of each community                  #
  # membership = an nx1 vector indicating the community membership for each node                  #
  #                                                                                               #
  # Outputs:                                                                                      #
  #          I = the observed information matrix which has size (B(B+1)/2 + n)x(B(B+1)/2 + n)     #
  #################################################################################################
  # get data related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  delta.mat <- (Adj == 0) * 1 # delta matrix indicating which elements of the adjacency matrix are 0
  
  # get the parameters 
  # little matrix for lambda (used for blockwise iterations)
  lambda.mat <- matrix(0, nrow = B, ncol = B)
  lambda.mat[lower.tri(lambda.mat, diag = TRUE)] <- param[1:num.blocks]
  lambda.mat <- lambda.mat + t(lambda.mat)
  diag(lambda.mat) <- diag(lambda.mat)/2
  num.na <- sum(is.na(lambda.mat))/2
  
  # big matrix for lambda (used for cellwise iterations)
  lambda.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      lambda.mat.BIG[i,j] <- lambda.mat[indx.r,indx.c]
    }
  }
  
  # vector for thetas
  theta.vec <- param[(num.blocks+1):length(param)]
  
  # matrix for thetas
  theta.mat <- matrix(theta.vec, ncol = 1) %*% matrix(theta.vec, nrow = 1)
  
  # Lambda Block
  lambda.deriv.vec <- rep(0, num.blocks-num.na)
  block.count <- 1
  for(i in 1:B){
    for(j in i:B){
      indx.i = which(membership == i)
      indx.j = which(membership == j)
      # get block-specific values of lambda
      lambda <- lambda.mat[i,j]
      if(!is.na(lambda)){# only do something if lambda is not NA
        # derivatives for lambda
        a <- Adj[indx.i, indx.j]
        delta <- delta.mat[indx.i, indx.j]
        theta1theta2 <- theta.mat[indx.i, indx.j]
        lambda.deriv <- (1 - delta) * (-a / lambda ^ 2 - theta1theta2 ^ 2 * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + theta1theta2 ^ 2 * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2)
        if(i == j){ # diagonal block gradients
          lambda.deriv.vec[block.count] <- sum(lambda.deriv[upper.tri(lambda.deriv)])  
        }else{ # off-diagonal block gradients
          lambda.deriv.vec[block.count] <- sum(lambda.deriv)  
        }
        block.count <- block.count + 1
      }
    }
  }
  lambda.block <- diag(lambda.deriv.vec)
  
  # Lamba-Theta Block
  lambdatheta.block <- matrix(0, nrow = num.blocks-num.na, ncol = n)
  block.count <- 1
  for(i in 1:B){
    for(j in i:B){
      indx.i <- which(membership == i)
      indx.j <- which(membership == j)
      lambda <- lambda.mat[i,j]
      if(!is.na(lambda)){# only do something if lambda is not NA
        delta.temp <- delta.mat[indx.i, indx.j]
        theta1theta2.temp <- theta.mat[indx.i, indx.j]
        if(i == j){
          for(k in indx.i){
            theta_k <- theta.vec[k]
            delta <- matrix(0, nrow = length(indx.i), ncol = length(indx.i)) #delta T
            delta[which(indx.i == k),] <- delta.temp[which(indx.i == k),]
            delta[,which(indx.i == k)] <- delta.temp[,which(indx.i == k)]
            theta1theta2 <- matrix(0, nrow = length(indx.i), ncol = length(indx.i)) #theta1theta2 T
            theta1theta2[which(indx.i == k),] <- theta1theta2.temp[which(indx.i == k),]
            theta1theta2[,which(indx.i == k)] <- theta1theta2.temp[,which(indx.i == k)]
            lambdatheta.deriv <- (1 - delta) * (-(theta1theta2/theta_k) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) - (theta1theta2^2 / theta_k) * lambda * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + (theta1theta2^2 / theta_k) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2 * lambda)
            lambdatheta.deriv[is.nan(lambdatheta.deriv)] <- 0
            lambdatheta.block[block.count,k] <- 0.5*sum(lambdatheta.deriv)
          }
        }else{
          for(k in indx.i){ #rows
            theta_k <- theta.vec[k]
            delta <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #delta Horiz  
            delta[which(indx.i == k),] <- delta.temp[which(indx.i == k),]
            theta1theta2 <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #theta1theta2 Horiz
            theta1theta2[which(indx.i == k),] <- theta1theta2.temp[which(indx.i == k),]
            lambda.deriv <- (1 - delta) * (-(theta1theta2/theta_k) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) - (theta1theta2^2 / theta_k) * lambda * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + (theta1theta2^2 / theta_k) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2 * lambda)
            lambdatheta.deriv[is.nan(lambdatheta.deriv)] <- 0
            lambdatheta.block[block.count,k] <- sum(lambdatheta.deriv)
          }
          for(k in indx.j){ #columns
            theta_k <- theta.vec[k]
            delta <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #delta Vert
            delta[,which(indx.i == k)] <- delta.temp[,which(indx.i == k)]
            theta1theta2 <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #theta1theta2 Vert
            theta1theta2[,which(indx.i == k)] <- theta1theta2.temp[,which(indx.i == k)]
            lambda.deriv <- (1 - delta) * (-(theta1theta2/theta_k) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) - (theta1theta2^2 / theta_k) * lambda * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + (theta1theta2^2 / theta_k) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2 * lambda)
            lambdatheta.deriv[is.nan(lambdatheta.deriv)] <- 0
            lambdatheta.block[block.count,k] <- sum(lambdatheta.deriv)
          }
        }
        block.count <- block.count + 1
      }
    }
  }
  
  # Theta Block Off-Diagonal
  #lambda.mat.BIG[is.na(lambda.mat.BIG)] <- 0 # replaces lambda NA's with 0's
  theta.block <- (1 - delta.mat) * (-lambda.mat.BIG * exp(lambda * theta.mat) / (exp(lambda.mat.BIG * theta.mat) - 1) - lambda.mat.BIG ^ 2 * theta.mat * exp(lambda.mat.BIG * theta.mat) / (exp(lambda.mat.BIG * theta.mat) - 1) + lambda.mat.BIG ^ 2 * theta.mat * exp(lambda.mat.BIG * theta.mat) ^ 2 / (exp(lambda.mat.BIG * theta.mat) - 1) ^ 2)
  theta.block[is.na(theta.block)] <- 0
  
  # Theta Block Diagonal
  for(k in 1:n){
    theta_k <- theta.vec[k]
    # Adjacency T
    a <- matrix(0, nrow = n, ncol = n)
    a[k,] <- Adj[k,]
    a[,k] <- Adj[,k]
    # Delta T
    delta <- matrix(0, nrow = n, ncol = n)
    delta[k,] <- delta.mat[k,]
    delta[,k] <- delta.mat[,k]
    # Lambda T
    lambda <- matrix(0, nrow = n, ncol = n)
    lambda[k,] <- lambda.mat.BIG[k,]
    lambda[,k] <- lambda.mat.BIG[,k]
    # theta1theta2 T
    theta1theta2 <- matrix(0, nrow = n, ncol = n)
    theta1theta2[k,] <- theta.mat[k,]
    theta1theta2[,k] <- theta.mat[,k]
    # derivatives with respect to theta_k
    theta.deriv <- (1 - delta) * (-a / theta_k ^ 2 - lambda ^ 2 * (theta1theta2^2 / theta_k^2) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + lambda ^ 2 * (theta1theta2^2 / theta_k^2) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2)
    theta.deriv[is.nan(theta.deriv)] <- 0
    theta.block[k,k] <- 0.5*sum(theta.deriv, na.rm = TRUE) #divide by 2 here because we only need to sum the upper triangle
  }
  
  I <- -rbind(cbind(lambda.block, lambdatheta.block), cbind(t(lambdatheta.block), theta.block))
  
  return(I)
}
# </END Observed Information Function>


# <BEGIN Expected Information Function>
exp.infomat <- function(param, exp.Adj, n.comm, membership){
  #################################################################################################
  #   This function computes the observed information matrix with the hurdle DCSBM. Note that     #
  #   this ignores the pi's since they can be estimated separately.                               #
  #                                                                                               #
  # Inputs:                                                                                       #
  #      param = a (B(B+1) + n)x1 parameter vector (the first B(B+1)/2 entries are devoted to the #
  #              pi'sm the second B(B+1)/2 entries are devoted to the lambda's and the final n    #
  #              entries are devoted to the theta's)                                              #
  #    exp.Adj = the nxn expected adjacency matrix                                                #
  #     n.comm = a Bx1 vector whose elements indicate the size of each community                  #
  # membership = an nx1 vector indicating the community membership for each node                  #
  #                                                                                               #
  # Outputs:                                                                                      #
  #          J = the expected information matrix which has size (B(B+1)/2 + n)x(B(B+1)/2 + n)     #
  #################################################################################################
  # get data related quantities
  B <- length(n.comm) # number of communities
  n <- sum(n.comm) # numbers of nodes
  num.blocks <- B*(B+1)/2 # number of unique community "blocks"
  
  # get the parameters 
  # little pi matrix (needed for the big pi matrix)
  pi.mat <- matrix(0, nrow = B, ncol = B)
  pi.mat[lower.tri(pi.mat, diag = TRUE)] <- param[1:num.blocks]
  pi.mat <- pi.mat + t(pi.mat)
  diag(pi.mat) <- diag(pi.mat)/2
  # big pi matrix 
  pi.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      pi.mat.BIG[i,j] <- pi.mat[indx.r,indx.c]
    }
  }
  # little matrix for lambda (used for blockwise iterations)
  lambda.mat <- matrix(0, nrow = B, ncol = B)
  lambda.mat[lower.tri(lambda.mat, diag = TRUE)] <- param[(num.blocks+1):(2*num.blocks)]
  lambda.mat <- lambda.mat + t(lambda.mat)
  diag(lambda.mat) <- diag(lambda.mat)/2
  num.na <- sum(is.na(lambda.mat))/2
  
  # big matrix for lambda (used for cellwise iterations)
  lambda.mat.BIG <- matrix(0, nrow = n, ncol = n)
  for(i in 1:n){
    for(j in 1:n){
      indx.r <- membership[i]
      indx.c <- membership[j]
      lambda.mat.BIG[i,j] <- lambda.mat[indx.r,indx.c]
    }
  }
  
  # vector for thetas
  theta.vec <- param[(2*num.blocks+1):length(param)]
  
  # matrix for thetas
  theta.mat <- matrix(theta.vec, ncol = 1) %*% matrix(theta.vec, nrow = 1)
  
  # Lambda Block
  lambda.deriv.vec <- rep(0, num.blocks-num.na)
  block.count <- 1
  for(i in 1:B){
    for(j in i:B){
      indx.i = which(membership == i)
      indx.j = which(membership == j)
      # get block-specific values of lambda
      lambda <- lambda.mat[i,j]
      if(!is.na(lambda)){# only do something if lambda is not NA
        # derivatives for lambda
        a <- exp.Adj[indx.i, indx.j]
        pi <- pi.mat.BIG[indx.i, indx.j]
        theta1theta2 <- theta.mat[indx.i, indx.j]
        lambda.deriv <- (1 - pi) * (-a / lambda ^ 2 - theta1theta2 ^ 2 * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + theta1theta2 ^ 2 * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2)
        if(i == j){ # diagonal block gradients
          lambda.deriv.vec[block.count] <- sum(lambda.deriv[upper.tri(lambda.deriv)])  
        }else{ # off-diagonal block gradients
          lambda.deriv.vec[block.count] <- sum(lambda.deriv)  
        }
        block.count <- block.count + 1
      }
    }
  }
  lambda.block <- diag(lambda.deriv.vec)
  
  # Lamba-Theta Block
  lambdatheta.block <- matrix(0, nrow = num.blocks-num.na, ncol = n)
  block.count <- 1
  for(i in 1:B){
    for(j in i:B){
      indx.i <- which(membership == i)
      indx.j <- which(membership == j)
      lambda <- lambda.mat[i,j]
      if(!is.na(lambda)){# only do something if lambda is not NA
        pi.temp <- pi.mat.BIG[indx.i, indx.j]
        theta1theta2.temp <- theta.mat[indx.i, indx.j]
        if(i == j){
          for(k in indx.i){
            theta_k <- theta.vec[k]
            pi <- matrix(0, nrow = length(indx.i), ncol = length(indx.i)) #pi T
            pi[which(indx.i == k),] <- pi.temp[which(indx.i == k),]
            pi[,which(indx.i == k)] <- pi.temp[,which(indx.i == k)]
            theta1theta2 <- matrix(0, nrow = length(indx.i), ncol = length(indx.i)) #theta1theta2 T
            theta1theta2[which(indx.i == k),] <- theta1theta2.temp[which(indx.i == k),]
            theta1theta2[,which(indx.i == k)] <- theta1theta2.temp[,which(indx.i == k)]
            lambdatheta.deriv <- (1 - pi) * (-(theta1theta2/theta_k) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) - (theta1theta2^2 / theta_k) * lambda * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + (theta1theta2^2 / theta_k) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2 * lambda)
            lambdatheta.deriv[is.nan(lambdatheta.deriv)] <- 0
            lambdatheta.block[block.count,k] <- 0.5*sum(lambdatheta.deriv)
          }
        }else{
          for(k in indx.i){ #rows
            theta_k <- theta.vec[k]
            pi <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #pi Horiz  
            pi[which(indx.i == k),] <- pi.temp[which(indx.i == k),]
            theta1theta2 <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #theta1theta2 Horiz
            theta1theta2[which(indx.i == k),] <- theta1theta2.temp[which(indx.i == k),]
            lambda.deriv <- (1 - pi) * (-(theta1theta2/theta_k) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) - (theta1theta2^2 / theta_k) * lambda * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + (theta1theta2^2 / theta_k) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2 * lambda)
            lambdatheta.deriv[is.nan(lambdatheta.deriv)] <- 0
            lambdatheta.block[block.count,k] <- sum(lambdatheta.deriv)
          }
          for(k in indx.j){ #columns
            theta_k <- theta.vec[k]
            pi <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #pi Vert
            pi[,which(indx.i == k)] <- pi.temp[,which(indx.i == k)]
            theta1theta2 <- matrix(0, nrow = length(indx.i), ncol = length(indx.j)) #theta1theta2 Vert
            theta1theta2[,which(indx.i == k)] <- theta1theta2.temp[,which(indx.i == k)]
            lambda.deriv <- (1 - pi) * (-(theta1theta2/theta_k) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) - (theta1theta2^2 / theta_k) * lambda * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + (theta1theta2^2 / theta_k) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2 * lambda)
            lambdatheta.deriv[is.nan(lambdatheta.deriv)] <- 0
            lambdatheta.block[block.count,k] <- sum(lambdatheta.deriv)
          }
        }
        block.count <- block.count + 1
      }
    }
  }
  
  # Theta Block Off-Diagonal
  #lambda.mat.BIG[is.na(lambda.mat.BIG)] <- 0 # replaces lambda NA's with 0's
  theta.block <- (1 - pi.mat.BIG) * (-lambda.mat.BIG * exp(lambda * theta.mat) / (exp(lambda.mat.BIG * theta.mat) - 1) - lambda.mat.BIG ^ 2 * theta.mat * exp(lambda.mat.BIG * theta.mat) / (exp(lambda.mat.BIG * theta.mat) - 1) + lambda.mat.BIG ^ 2 * theta.mat * exp(lambda.mat.BIG * theta.mat) ^ 2 / (exp(lambda.mat.BIG * theta.mat) - 1) ^ 2)
  theta.block[is.na(theta.block)] <- 0
  
  # Theta Block Diagonal
  for(k in 1:n){
    theta_k <- theta.vec[k]
    # Adjacency T
    a <- matrix(0, nrow = n, ncol = n)
    a[k,] <- exp.Adj[k,]
    a[,k] <- exp.Adj[,k]
    # Pi T
    pi <- matrix(0, nrow = n, ncol = n)
    pi[k,] <- pi.mat.BIG[k,]
    pi[,k] <- pi.mat.BIG[,k]
    # Lambda T
    lambda <- matrix(0, nrow = n, ncol = n)
    lambda[k,] <- lambda.mat.BIG[k,]
    lambda[,k] <- lambda.mat.BIG[,k]
    # theta1theta2 T
    theta1theta2 <- matrix(0, nrow = n, ncol = n)
    theta1theta2[k,] <- theta.mat[k,]
    theta1theta2[,k] <- theta.mat[,k]
    # derivatives with respect to theta_k
    theta.deriv <- (1 - pi) * (-a / theta_k ^ 2 - lambda ^ 2 * (theta1theta2^2 / theta_k^2) * exp(lambda * theta1theta2) / (exp(lambda * theta1theta2) - 1) + lambda ^ 2 * (theta1theta2^2 / theta_k^2) * exp(lambda * theta1theta2) ^ 2 / (exp(lambda * theta1theta2) - 1) ^ 2)
    theta.deriv[is.nan(theta.deriv)] <- 0
    theta.block[k,k] <- 0.5*sum(theta.deriv, na.rm = TRUE) #divide by 2 here because we only need to sum the upper triangle
  }
  
  J <- -rbind(cbind(lambda.block, lambdatheta.block), cbind(t(lambdatheta.block), theta.block))
  
  return(J)
}
# </END Expected Information Function>


# <BEGIN NA Inserter>
insert.na <- function(n, x, idx){
  # Inputs:
  #    n = length of original vector from which elements were deleted
  #    x = the new vector into which you wish to insert NAs
  #  idx = a vectors of indices (relative to the original vector) where you want NAs. These indices must be sorted from smallest to largest.
  for(i in 1:length(idx)){
    x <- append(x = x, values = NA, after = length(x) - (n - idx[length(idx) + 1 - i]))
  }
  return(x)
}
# </END NA Inserter>