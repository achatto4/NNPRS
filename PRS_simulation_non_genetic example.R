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

# Load necessary libraries
library(Matrix)
library(Rcpp)
library(ggplot2)

# Function to generate synthetic LD matrix and GWAS summary statistics
generate_synthetic_data <- function(n, p, beta) {
  X <- matrix(rnorm(n * p), n, p)
  R <- crossprod(X) / n
  y <- X %*% beta
  r <- crossprod(X, y) / n
  return(list(X = X, R = R, r = r, y = y))
}

# Parameters
p <- 1000  # Number of SNPs
n_EUR <- 1000  # Sample size for EUR
n_SAS <- 200  # Sample size for SAS
n_EAS <- 200  # Sample size for EAS
n_AFR <- 100  # Sample size for AFR (main population)

# Set parameters
percentage_nonzero <- 0.001  # 1% non-zero SNPs
num_nonzero <- ceiling(p * percentage_nonzero)

# Initialize beta with zeros
beta1 <- rep(0, p)

# Generate non-zero indices
set.seed(123)
nonzero_indices <- sample(1:p, num_nonzero)

# Assign random values to non-zero indices
beta1[nonzero_indices] <- rnorm(num_nonzero, mean = 0, sd = 0.001)

# Generate synthetic data for auxiliary populations
data_EUR <- generate_synthetic_data(n_EUR, p, beta1)
data_SAS <- generate_synthetic_data(n_SAS, p, beta1)
data_EAS <- generate_synthetic_data(n_EAS, p, beta1)

# Generate synthetic data for main population with original beta
data_AFR <- generate_synthetic_data(n_AFR, p, beta1)

# Introduce jitter to beta for AFR population
set.seed(1)
beta_AFR <- beta1 + rnorm(p, mean = 0, sd = 0.0001)

data_AFR_jittered <- generate_synthetic_data(n_AFR, p, beta_AFR)

# Extract necessary matrices
R_AFR <- data_AFR$R
r_AFR <- data_AFR$r
R_AFR_jittered <- data_AFR_jittered$R
r_AFR_jittered <- data_AFR_jittered$r

# Combine data into lists
R_list <- list(data_EUR$R, data_SAS$R, data_EAS$R)
r_list <- list(data_EUR$r, data_SAS$r, data_EAS$r)
n_list <- c(n_EUR, n_SAS, n_EAS)

alpha = 0.1
eta = 1
# Algorithm parameters
alpha1 <-alpha
alpha2 <-alpha
alpha3 <-alpha
alpha4 <-alpha
eta_l <-eta
eta_m <-eta
max_iter <- 1000

# Load Rcpp function
library(Rcpp)
sourceCpp("grad_func.cpp")

# Run gradient descent algorithm

# res<- gradient_descent_transfer_learning(
#     n_AFR, r_AFR, R_AFR, n_list, r_list, R_list,
#     alpha1, alpha2, alpha3, alpha4, eta_l, eta_m, max_iter
#   )

# res_rcpp <- gradient_descent_transfer_learning_rcpp(
#   n_AFR, r_AFR, R_AFR, n_list, r_list, R_list,
#   alpha1, alpha2, alpha3, alpha4, eta_l, eta_m, max_iter
# )

res_rcpp_jittered <- gradient_descent_transfer_learning_rcpp(
  n_AFR, r_AFR_jittered, R_AFR_jittered, n_list, r_list, R_list,
  alpha1, alpha2, alpha3, alpha4, eta_l, eta_m, max_iter
)

res_rcpp_P4 <- gradient_descent_main_P4(
  n_AFR, r_AFR_jittered, R_AFR_jittered,
  alpha1, alpha2, eta_m, max_iter
)

res_rcpp_single <- gradient_descent_main_only(
  n_AFR, r_AFR_jittered, R_AFR_jittered,
  alpha1, alpha2, eta_m, max_iter
)

#Compute PRS
#Compute the 99th percentile threshold of absolute values
#threshold <- quantile(abs(res_rcpp_jittered$hat_beta), 0)
#threshold <- quantile(abs(res$hat_beta), 0.99)
# Set values below the threshold to 0
#res_rcpp_jittered$hat_beta[abs(res_rcpp_jittered$hat_beta) < threshold] <- 0
# res$hat_beta[abs(res$hat_beta) < threshold] <- 0

PRS_original <- data_AFR_jittered$X %*% beta_AFR
#PRS_jittered <- data_AFR_jittered$X %*% res_rcpp_jittered$hat_beta
PRS_jittered <- data_AFR_jittered$X %*% res_rcpp_jittered$hat_beta
PRS_single <- data_AFR_jittered$X %*% res_rcpp_single$hat_beta
PRS_sinP4 <- data_AFR_jittered$X %*% res_rcpp_P4$hat_beta
rank(PRS_original); rank(PRS_jittered); rank(PRS_single); rank(PRS_sinP4)
kendall_tau <- cor(PRS_original, PRS_jittered, method = "kendall")
print(kendall_tau)

# Fit linear models
model_original <- lm(data_AFR_jittered$y ~ PRS_original)
model_jittered <- lm(data_AFR_jittered$y ~ PRS_jittered)
model_single <- lm(data_AFR_jittered$y ~ PRS_single)
model_sinP4 <- lm(data_AFR_jittered$y ~ PRS_sinP4)
# Compute R^2 values
r2_original <- summary(model_original)$r.squared
r2_jittered <- summary(model_jittered)$r.squared
r2_single <- summary(model_single)$r.squared
r2_sinP4 <- summary(model_sinP4)$r.squared
# Print results
cat("R^2 for y ~ PRS_original:", r2_original, "\n")
cat("R^2 for y ~ PRS_jittered:", r2_jittered, "\n")
cat("R^2 for y ~ PRS_single:", r2_single, "\n")
cat("R^2 for y ~ PRS_sinP4:", r2_sinP4, "\n")
# Plot PRS distributions
ggplot() +
  geom_density(aes(PRS_original), fill = "red", alpha = 0.5) +
  geom_density(aes(PRS_jittered), fill = "blue", alpha = 0.5) +
  ggtitle("PRS Distribution: Original vs Jittered") +
  xlab("PRS Score") +
  theme_minimal()

ggplot() +
  geom_density(aes(PRS_original), fill = "red", alpha = 0.5) +
  geom_density(aes(PRS_single), fill = "blue", alpha = 0.5) +
  ggtitle("PRS Distribution: Original vs single") +
  xlab("PRS Score") +
  theme_minimal()

ggplot() +
  geom_density(aes(PRS_original), fill = "red", alpha = 0.5) +
  geom_density(aes(PRS_sinP4), fill = "blue", alpha = 0.5) +
  ggtitle("PRS Distribution: Original vs sinP4") +
  xlab("PRS Score") +
  theme_minimal()

beta_AFR[nonzero_indices]
res_rcpp_jittered$hat_beta[nonzero_indices]
res_rcpp_single$hat_beta[nonzero_indices]
res_rcpp_P4$hat_beta[nonzero_indices]
