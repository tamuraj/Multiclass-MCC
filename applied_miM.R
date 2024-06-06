###############################################################
# Point estimates and asymptotic confidence intervals for miM #
###############################################################

# The code derives point estimates and asymptotic confidence intervals for miM. The function arguments are shown below.

# prob_matri: Matrix of probabilities
# n         : samplesize
# conf.level: confidence interval level


compute_microMcc_CI <- function(prob_matrix, n, conf.level) {
  
  
  micro_averaged_mcc <- function(prob_matrix) {
    r <- nrow(prob_matrix)
    di <- sum(diag(prob_matrix))
    micro_mcc <- (r * di - 1) / (r - 1)
    return(micro_mcc)
  }
  
  # Function to compute A_{ii}^{miM}
  compute_Aii_miM <- function(r) {
    Aii_miM <- r / (1 - r)
    return(Aii_miM)
  }
  
  # Function to compute A_{ij}^{miM}
  compute_Aij_miM <- function() {
    Aij_miM <- 0
    return(Aij_miM)
  }
  
  compute_gradient <- function(prob_matrix) {
    r <- nrow(prob_matrix)
    gradient <- numeric(r * r)
    
    for (i in 1:r) {
      for (j in 1:r) {
        if (i == j) {
          gradient[(i - 1) * r + j] <- compute_Aii_miM(r)
        } else {
          gradient[(i - 1) * r + j] <- compute_Aij_miM()
        }
      }
    }
    
    return(gradient)
  }
  
 
  mcc <- micro_averaged_mcc(prob_matrix)
  
  
  d <- compute_gradient(prob_matrix)
  pa <- c(t(prob_matrix))
  d1 <- c(d)
  p <- t(t(pa))
  dp <- diag(pa)
  ma <- rbind(d1)
  
  
  s1 <- ma %*% (dp - p %*% t(p)) %*% t(ma)
  
  
  z <- qnorm(1 - (1 - conf.level) / 2)
  v <- s1 / n 
  s <- sqrt(v)
  ci <- mcc + c(-1, 1) * z * s
  
  output <- c("miM" = mcc, "confidence_interval" = ci)
  return(output)
}

# Example 
prob_matrix <- matrix(c(
  0.1, 0.05, 0.05,
  0.1, 0.2, 0.1,
  0.05, 0.1, 0.15
), nrow = 3, byrow = TRUE)

prob_matrix <- prob_matrix / sum(prob_matrix)
compute_microMcc_CI(prob_matrix, 1000, 0.95)
