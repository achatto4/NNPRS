# Load required libraries
library(matrixStats)

# Define the gradient descent function for transfer learning with corrected scaling
gradient_descent_transfer_learning <- function(n0, r0, R0, nk_list, rk_list, Rk_list, alpha1, alpha2, alpha3, alpha4, eta_l, eta_m, max_iter) {
  # Initialize variables
  u_l <- alpha1 * rep(1, dim(R0)[1])
  v_l <- alpha2 * rep(1, dim(R0)[1])
  h_m <- alpha3 * rep(1, dim(R0)[1])
  g_m <- alpha4 * rep(1, dim(R0)[1])
  
  # Gradient descent for auxiliary data
  for (l in 0:max_iter) {
    grad_u <- rep(0, length(u_l))
    grad_v <- rep(0, length(v_l))
    
    for (k in 1:length(nk_list)) {
      nk <- nk_list[[k]]
      rk <- rk_list[[k]]
      Rk <- Rk_list[[k]]
      
      u_l_sq <- u_l * u_l
      v_l_sq <- v_l * v_l
      diff_sq <- u_l_sq - v_l_sq
      
      grad_u <- grad_u + (-4 * nk * rk * u_l + 4 * nk * (Rk %*% diff_sq) * u_l)
      grad_v <- grad_v + (4 * nk * rk * v_l - 4 * nk * (Rk %*% diff_sq) * v_l)
    }
    
    total_n <- sum(unlist(nk_list)) + n0
    u_l <- u_l - (eta_l / total_n) * grad_u
    v_l <- v_l - (eta_l / total_n) * grad_v
  }
  
  hat_u <- u_l
  hat_v <- v_l
  
  # Gradient descent for main data
  for (m in 0:max_iter) {
    u_hat_sq <- hat_u * hat_u
    v_hat_sq <- hat_v * hat_v
    grad_h <- rep(0, length(h_m))
    grad_g <- rep(0, length(g_m))
    
    h_m_sq <- h_m * h_m
    g_m_sq <- g_m * g_m
    diff_sq <- u_hat_sq - v_hat_sq + h_m_sq - g_m_sq
    
    grad_h <- (-4 * n0 * r0 * h_m + 4 * n0 * (R0 %*% diff_sq) * h_m)
    grad_g <- (4 * n0 * r0 * g_m - 4 * n0 * (R0 %*% diff_sq) * g_m)
    
    h_m <- h_m - (eta_m / n0) * grad_h
    g_m <- g_m - (eta_m / n0) * grad_g
  }
  
  hat_h <- h_m
  hat_g <- g_m
  
  # Define the output
  hat_beta <- hat_u * hat_u - hat_v * hat_v + hat_h * hat_h - hat_g * hat_g
  
  return(list(hat_u = hat_u, hat_v = hat_v, hat_h = hat_h, hat_g = hat_g, hat_beta = hat_beta))
}

# Load necessary library
library(Matrix)


# Function to generate synthetic LD matrix and GWAS summary statistics
generate_synthetic_data <- function(n, p, beta = beta1) {
  # Generate random LD matrix
  X <- matrix(rnorm(n * p), n, p)
  R <- crossprod(X) / n
  
  # Generate random GWAS summary statistics
  y <- X%*%beta
  r <- crossprod(X, y) / n
  
  return(list(R = R, r = r))
}

# Parameters
p <- 1000  # Number of SNPs (keeping it small for simplicity)
n_EUR <- 10000  # Sample size for EUR
n_SAS <- 2000  # Sample size for SAS
n_EAS <- 2000  # Sample size for EAS
n_AFR <- 100  # Sample size for AFR (main population)

# Set parameters
percentage_nonzero <- 0.01 # 1% non-zero
num_nonzero <- ceiling(p * percentage_nonzero) # Calculate the number of non-zero elements

# Initialize beta with zeros
beta1 <- rep(0, p)

# Generate non-zero indices
set.seed(123) # For reproducibility
nonzero_indices <- sample(1:p, num_nonzero)

# Assign random values to non-zero indices
beta1[nonzero_indices] <- rnorm(num_nonzero, mean = 0, sd = 0.5)


# Generate synthetic data for auxiliary populations
data_EUR <- generate_synthetic_data(n_EUR, p)
data_SAS <- generate_synthetic_data(n_SAS, p)
data_EAS <- generate_synthetic_data(n_EAS, p)

# Generate synthetic data for main population
data_AFR <- generate_synthetic_data(n_AFR, p)

# Extract R and r for each population
R_EUR <- data_EUR$R
r_EUR <- data_EUR$r

R_SAS <- data_SAS$R
r_SAS <- data_SAS$r

R_EAS <- data_EAS$R
r_EAS <- data_EAS$r

R_AFR <- data_AFR$R
r_AFR <- data_AFR$r

# Combine data into lists for easier processing
R_list <- list(R_EUR, R_SAS, R_EAS)
r_list <- list(r_EUR, r_SAS, r_EAS)
n_list <- c(n_EUR, n_SAS, n_EAS)

# Set algorithm parameters
alpha1 <- 0.01
alpha2 <- 0.01
alpha3 <- 0.01
alpha4 <- 0.01
eta_l <- 0.01
eta_m <- 0.01
max_iter <- 10000

# Run the gradient descent algorithm
result <- gradient_descent_transfer_learning(n_AFR, r_AFR, R_AFR, n_list, r_list, R_list, alpha1, alpha2, alpha3, alpha4, eta_l, eta_m, max_iter)
# Load Rcpp
library(Rcpp)

sourceCpp("grad_func.cpp")
res_rcpp <- gradient_descent_transfer_learning_rcpp(
    n_AFR, r_AFR, R_AFR, n_list, r_list, R_list, 
    alpha1, alpha2, alpha3, alpha4, eta_l, eta_m, max_iter
  )
# Print the results
cat("hat_u:\n")
print(result$hat_u)
cat("hat_v:\n")
print(result$hat_v)
cat("hat_h:\n")
print(result$hat_h)
cat("hat_g:\n")
print(result$hat_g)
cat("hat_beta:\n")
print(result$hat_beta)

# Check gradients convergence
cat("Final gradient for h_m:\n")
print(result$hat_h - h_m)
cat("Final gradient for g_m:\n")
print(result$hat_g - g_m)

# Beta values for AFR
cat("Beta values for AFR:\n")
print(result$hat_beta)
cat("true Beta values for AFR:\n")
beta1

res_rcpp$hat_beta
cat("Beta values for AFR just using AFR:\n")
solve(R_AFR,r_AFR)