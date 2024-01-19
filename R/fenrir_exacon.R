################################################################################

# fenrir -- exacon(), a function to replicate EXACON

################################################################################

# exacon -----------------------------------------------------------------------

# This function is designed to replicate the behavior of the EXACON program
# (Bergman & El-Khouri, 1987, https://doi.org/10.1177/0013164487471024).

# To do:
# - Test to see what breaks the function
# - Allow for 0 frequencies

exacon <- function(f) {
  
  # Extract values
  
  r <- rowSums(f)
  
  c <- colSums(f)
  
  sum_grid <- expand.grid(r, c)
  
  N <- sum(f)
  
  e <- matrix((sum_grid[, 1] * sum_grid[, 2]) / N, nrow = nrow(f))
  
  p <- matrix((sum_grid / N)[, 1] * (sum_grid / N)[, 2], nrow = nrow(f))
  
  q <- 1 - p
  
  # Internal functions
  
  ## Calculate the p-value for the binomial test (complete independence model)
  
  p_binomial <- function(f, e, i, N, p, q) {
    
    cfn <- choose(N, f[i])
    
    p_f <- cfn * p[i]^f[i] * q[i]^(N - f[i])
    
    if (f[i] > e[i]) {
      
      p_i    <- rep(NA, length(f[i]:N) + 1)
      
      p_i[1] <- p_f
      
      for (j in 2:(length(f[i]:N) + 1)) {
        
        x      <- f[i] + j - 2
        
        g <- (p[i] * (N - x)) / (q[i] * (x + 1))
        
        p_i[j] <- p_i[j - 1] * g
        
      }
      
      sum(p_i)
      
    } else {
      
      p_i    <- rep(NA, length(0:f[i]))
      
      p_i[1] <- p_f
      
      for (j in 1:f[i]) {
        
        x      <- f[i] - j + 1
        
        g <- (x * q[i]) / (p[i] * (N - x + 1))
        
        p_i[j + 1] <- p_i[j] * g
        
      }
      
      sum(p_i)
      
    }
    
  }
  
  ## Calculate p-value for the hypergeometric test (fixed margins model)
  
  p_hyper <- function(f, e, i, N, p, q, sum_grid) {
    
    r <- sum_grid[i, 1]
    c <- sum_grid[i, 2]
    
    crn     <- choose(N, r)
    crn_sub <- choose(N - c, r - f[i])
    cfc     <- choose(c, f[i])
    
    p_f <- (cfc * crn_sub) / crn
    
    if (f[i] > e[i]) {
      
      p_i    <- rep(NA, length(f[i]:N) + 1)
      
      p_i[1] <- p_f
      
      for (j in 2:(length(f[i]:N) + 1)) {
        
        x      <- f[i] + j - 2
        
        g <- ((c - x) * (r - x)) / ((x + 1) * (N - c - r + x + 1))
        
        p_i[j] <- p_i[j - 1] * g
        
      }
      
      sum(p_i)
      
    } else {
      
      p_i    <- rep(NA, length(0:f[i]))
      
      p_i[1] <- p_f
      
      for (j in 1:f[i]) {
        
        x      <- f[i] - j + 1
        
        g <- (x * (N - c - r + x)) / ((c - x + 1) * (r - x + 1))
        
        p_i[j + 1] <- p_i[j] * g
        
      }
      
      sum(p_i)
      
    }
    
  }
  
  # Calculate binomial and hypergeometric p-values
  
  index <- 1:length(f)
  
  pb <- matrix(rep(NA, length(f)), nrow = nrow(f))
  ph <- matrix(rep(NA, length(f)), nrow = nrow(f))
  
  for (i in index) {
    
    pb[i] <- p_binomial(f, e, index[i], N, p, q)
    ph[i] <- p_hyper(f, e, index[i], N, p, q, sum_grid)
    
  }

  # Create output
  
  out <- list(
    "observed" = f,
    "expected" = e,
    "bin_prob" = pb,
    "hyp_prob" = ph
  )
  
  return(out)
  
}
