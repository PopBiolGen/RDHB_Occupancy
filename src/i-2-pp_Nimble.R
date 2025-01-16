# A nimble model to estimate parameters of the pp model.

library(nimble)

# simulate the data
source("src/i-1-simulate-pp.R")

# Nimble function definitions
pairwise_distance_nf <- nimbleFunction(
  run = function(A = double(2), B = double(2)) {
    # Input: A and B are matrices (2D arrays)
    # Output: Pairwise distance matrix
    
    # Declare return type
    returnType(double(2))
    
    # Get dimensions of A and B
    nA <- dim(A)[1]  # Number of rows in A
    nB <- dim(B)[1]  # Number of rows in B
    
    # Initialize the distance matrix
    D <- matrix(0, nrow = nA, ncol = nB)
    
    # Compute the pairwise distances
    for (i in 1:nA) {
      for (j in 1:nB) {
        # d_ij^2 = ||a_i||^2 + ||b_j||^2 - 2 * dot(a_i, b_j)
        D[i, j] <- pow((A[i, 1]-B[j, 1])^2 + (A[i, 2]-B[j, 2])^2, 1/2)
      }
    }
    
    return(D)  # Return the pairwise distance matrix
  }
)

kernel_prod_nf <- nimbleFunction(
  run = function(distmat = double(2), z = double(1), k = double(0)) {
    # Input: A and B are distance matrix and scalar for k
    # Output: Probability matrix
    
    # Declare return type
    returnType(double(1))
    
    # Get dimensions of A
    dim.dm <- dim(distmat)  # Dimensions of distmat

    # Initialize the probability matrix
    P <- numeric(dim.dm[1], 1)
    
    # Compute the probability vector (sum over the transformed distance matrix)
    for (i in 1:dim.dm[1]) {
      for (j in 1:dim.dm[2]) {
        # probability of not present from any source site
        P[i] <- P[i] * (z[j] * (1-exp(-distmat[i,j]/k)) + 1-z[j]) # product over (real) colonies
      }
    }
    P[1:dim.dm[1]] <- 1 - P[1:dim.dm[1]] # probability present
    
    return(P)  # Return the probability matrix
  }
)

pp.code <- nimbleCode(
 {
  #priors
  psi ~ dunif(0, 1) # probability of a DA colony being real
  r0 ~ dunif(500, x.max/2) # extent of the invasion
  alpha.det ~ dnorm(0, 1/4) # intercept for log-odds of detection
  beta.det ~ dnorm(0, 1/4) # slope of some linear effect on detection
  k ~ dunif(100, 10000) # parameter affecting spread of the density kernel
  
  # place ground 0
  g0[1, 1] ~ dunif(g0.x.min, g0.x.max) # ground zero x
  g0[1, 2] ~ dunif(g0.y.min, g0.y.max) # ground zero y
  g0[2, 1] <- 0 # these just to stop Nimble dropping dimensions
  g0[2, 2] <- 0
  
  # place colonies
  for (jj in 1:JJ){
    c.j0[jj, 1] ~ dunif(x.min, x.max) # place initial colonies into X Y matrix
    c.j0[jj, 2] ~ dunif(y.min, y.max)
    
    # which colonies are real?
    da[jj] ~ dbern(psi) # fraction of data augmented colonies
  }
  
  # find colonies within extent radius
  # Compute the pairwise distances
  in.ex[1:JJ] <- step(pairwise_distance_nf(c.j0[,], g0[,])[1:JJ,1] - r0)
  z0[1:JJ] <- da[1:JJ]*in.ex[1:JJ] # final z-variable defining which DA colonies are 'real'
  
  # calculate pairwise distance between sites and colonies
  pd[1:II, 1:JJ] <- pairwise_distance_nf(s.0[1:II, 1:2], c.j0[1:JJ, 1:2])
  # convert to a probability matrix based on the kernel
  prob.ij[1:II] <- kernel_prod_nf(pd[1:II, 1:JJ], k)
  
  # survey-level detection
  logit(det[1:II]) <- alpha.det + beta.det*sur.lev.var[1:II] # linear model on detection
  
  # calculate likelihoods
  for (ii in 1:II){
    obs.i[ii] ~ dbern(det[ii]*prob.ij[ii])
  }
  
 }
)

# Data
max.c <- ceiling(2*c.n) # maximum colonies (for data augmentation)

constant.list <- list(
  #nt = 1, # number of time steps
  x.min = 0,
  y.min = 0, # possible spatial extent of the species across all time, bounding box
  x.max = 20000,
  y.max = 20000, # in metres
  g0.x.min = 7500, # bounding box for possible location of ground zero
  g0.y.min = 7500,
  g0.x.max = 12500,
  g0.y.max = 12500,
  II = length(s.i0.x), # number of sites
  JJ = max.c # maximum number of colonies (data-augmentation approach)
)

data.list <- list(
  s.0 = s.0, # survey site locations (matrix of X and Y coordinates)
  sur.lev.var = sur.lev.var, # detection covariate (vector of length II)
  obs.i = obs.i
)

dim.list <- list(
  g0 = c(2, 2), # has to be 2x2 else Nimble makes it a vector
  c.j0 = c(max.c, 2),
  da = max.c
)

# initials
init.list <- list(g0 = matrix(c(8000, 0, 9000, 0), nrow = 2))

# parameters to monitor
params <- c("psi",
            "r0",
            "alpha.det",
            "beta.det",
            "k")

# mcmc settings
nb <- 5000
ni <- 10000
nc = 3

# the model
a.n <- nimbleModel(code = pp.code, 
                   inits = init.list, 
                   dimensions = dim.list,
                   constants = constant.list,
                   data = data.list)

a.n.c <- nimbleMCMC(code = pp.code, 
                  inits = init.list, 
                  dimensions = dim.list,
                  monitors = params, 
                  constants = constant.list,
                  data = data.list,
                  niter = ni, 
                  nburnin = nb, 
                  nchains = nc,
                  check = FALSE,
                  samplesAsCodaMCMC = TRUE)
