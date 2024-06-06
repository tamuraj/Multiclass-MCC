################################################################
# Point estimates and asymptotic confidence intervals for miM* #
################################################################

# The code derives point estimates and asymptotic confidence intervals for miM*. The function arguments are shown below.

# prob_matri: Matrix of probabilities
# n         : samplesize
# conf.level: confidence interval level

compute_Aii_miM_star <- function(prob_matrix, i, r) {
  p_i_dot <- sum(prob_matrix[i, ])
  p_dot_i <- sum(prob_matrix[, i])
  
  sum_pii <- sum(diag(prob_matrix))
  sum_pi_dot_p_dot_i <- sum(rowSums(prob_matrix) * colSums(prob_matrix))
  sum_p_dot_k_squared <- sum(colSums(prob_matrix)^2)
  sum_pi_dot_k_squared <- sum(rowSums(prob_matrix)^2)
  
  term1 <- (1 - p_i_dot - p_dot_i) / (sqrt(1 - sum_p_dot_k_squared) * sqrt(1 - sum_pi_dot_k_squared))
  
  term2_inner1 <- (sum_pii - sum_pi_dot_p_dot_i) * p_i_dot
  term2_inner2 <- (sum_pii - sum_pi_dot_p_dot_i) * p_dot_i
  term2_denom1 <- sqrt((1 - sum_p_dot_k_squared) * (1 - sum_pi_dot_k_squared)^3)
  term2_denom2 <- sqrt((1 - sum_pi_dot_k_squared) * (1 - sum_p_dot_k_squared)^3)
  
  term2 <- (term2_inner1 / term2_denom1) + (term2_inner2 / term2_denom2)
  
  Aii_miM_star <- term1 + term2
  
  return(Aii_miM_star)
}

compute_extendedMicroMcc_CI <- function(prob_matrix, n, conf.level) {
  
  # Function to calculate extended miM* 
  extended_micro_averaged_mcc <- function(prob_matrix) {
    r <- nrow(prob_matrix)
    
    sum_pii <- sum(diag(prob_matrix))
    sum_pi_dot_p_dot_i <- sum(rowSums(prob_matrix) * colSums(prob_matrix))
    
    numerator <- sum_pii - sum_pi_dot_p_dot_i
    denominator1 <- sqrt(1 - sum(colSums(prob_matrix)^2))
    denominator2 <- sqrt(1 - sum(rowSums(prob_matrix)^2))
    
    miM_star <- numerator / (denominator1 * denominator2)
    
    return(miM_star)
  }
  
  # Function to compute A_{ii}^{miM*}
  compute_Aii_miM_star <- function(prob_matrix, i, r) {
    p_i_dot <- sum(prob_matrix[i, ])
    p_dot_i <- sum(prob_matrix[, i])
    
    sum_pii <- sum(diag(prob_matrix))
    sum_pi_dot_p_dot_i <- sum(rowSums(prob_matrix) * colSums(prob_matrix))
    sum_p_dot_k_squared <- sum(colSums(prob_matrix)^2)
    sum_pi_dot_k_squared <- sum(rowSums(prob_matrix)^2)
    
    term1 <- (1 - p_i_dot - p_dot_i) / (sqrt(1 - sum_p_dot_k_squared) * sqrt(1 - sum_pi_dot_k_squared))
    
    term2_inner1 <- (sum_pii - sum_pi_dot_p_dot_i) * p_i_dot
    term2_inner2 <- (sum_pii - sum_pi_dot_p_dot_i) * p_dot_i
    term2_denom1 <- sqrt((1 - sum_p_dot_k_squared) * (1 - sum_pi_dot_k_squared)^3)
    term2_denom2 <- sqrt((1 - sum_pi_dot_k_squared) * (1 - sum_p_dot_k_squared)^3)
    
    term2 <- (term2_inner1 / term2_denom1) + (term2_inner2 / term2_denom2)
    
    Aii_miM_star <- term1 + term2
    
    return(Aii_miM_star)
  }
  
  # Function to compute A_{ij}^{miM*}
  compute_Aij_miM_star <- function(prob_matrix, i, j, r) {
    p_i_dot <- sum(prob_matrix[i,])
    p_dot_i <- sum(prob_matrix[,i])
    p_j_dot <- sum(prob_matrix[j,])
    p_dot_j <- sum(prob_matrix[,j])
    
    sum_pii <- sum(diag(prob_matrix))
    sum_p_dot_p_dot_i <- sum(rowSums(prob_matrix) * colSums(prob_matrix))
    sum_p_i_dot_k_squared <- sum(colSums(prob_matrix)^2)
    sum_p_dot_k_squared <- sum(rowSums(prob_matrix)^2)
    
    numerator1 <- -p_dot_i - p_j_dot
    denominator1 <- sqrt(1 - sum_p_dot_k_squared) * sqrt(1 - sum_p_i_dot_k_squared)
    
    numerator2 <- sum_pii - sum_p_dot_p_dot_i
    term1 <- p_i_dot / (sqrt(1 - sum_p_dot_k_squared)^3 * sqrt((1 - sum_p_i_dot_k_squared)))
    term2 <- p_dot_j / (sqrt(1 - sum_p_i_dot_k_squared)^3 * sqrt((1 - sum_p_dot_k_squared)))
    
    Aij_miM_star <- numerator1 / denominator1 + numerator2 * (term1 + term2)
    
    return(Aij_miM_star)
  }
  
  
  compute_gradient <- function(prob_matrix) {
    r <- nrow(prob_matrix)
    
    sum_pii <- sum(diag(prob_matrix))
    sum_pi_dot_p_dot_i <- sum(rowSums(prob_matrix) * colSums(prob_matrix))
    sum_p_dot_k_squared <- sum(colSums(prob_matrix)^2)
    sum_pi_dot_k_squared <- sum(rowSums(prob_matrix)^2)
    
    gradient <- numeric(r * r)
    
    for (i in 1:r) {
      for (j in 1:r) {
        pi_i_dot <- sum(prob_matrix[i, ])
        p_dot_i <- sum(prob_matrix[, i])
        
        if (i == j) {
          gradient[(i - 1) * r + j] <- compute_Aii_miM_star(
            prob_matrix, i, r
          )
        } else {
          gradient[(i - 1) * r + j] <- compute_Aij_miM_star(
            prob_matrix, i,j, r
          )
        }
      }
    }
    
    return(gradient)
  }
  
  # Calculate miM*
  miM_star <- extended_micro_averaged_mcc(prob_matrix)
  
 
  d <- compute_gradient(prob_matrix)
  pa <- c(t(prob_matrix))
  d1 <- c(d)
  p <- t(t(pa))
  dp <- diag(pa)
  ma <- rbind(d1)
  
  
  s1 <- ma %*% (dp - p %*% t(p)) %*% t(ma)
  
  # Calculate confidence interval
  z <- qnorm(1 - (1 - conf.level) / 2)
  v <- s1 / n 
  s <- sqrt(v)
  ci <- miM_star + c(-1, 1) * z * s
  
  output <- c("miM_star" = miM_star, "confidence_interval" = ci)
  return(output)
}
  


# Example 
prob_matrix <- matrix(c(
  0.1, 0.05, 0.05,
  0.1, 0.2, 0.1,
  0.05, 0.1, 0.15
), nrow = 3, byrow = TRUE)

compute_extendedMicroMcc_CI(prob_matrix, 1000, 0.95)



