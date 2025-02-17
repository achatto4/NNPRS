# Load required libraries
# Install 'bnlearn' package if not already installed for independence testing
if (!requireNamespace("bnlearn", quietly = TRUE)) {
  install.packages("bnlearn")
}
library(bnlearn)

load("data/simulated_genotype/s1/s1_9.RData")
#View(G_sim_fit)
#View(protein_network_sim)



# Function to identify root nodes
find_root_nodes <- function(data, instruments, level = 0.001) {
  # Inputs:
  #   data: A data frame containing nodes as columns and their respective values.
  #   instruments: A data frame with instrument variables corresponding to each node.
  
  n_nodes <- ncol(data)  # Number of nodes
  node_names <- colnames(data)
  root_set <- c()  # Initialize an empty set of roots
  
  # Step 2: Test independence for all node pairs
  for (i in 1:n_nodes) {
    is_root <- TRUE  # Assume node i is a root unless proven otherwise
    
    for (j in 1:n_nodes) {
      if (i != j) {
        # Step 2.1: Perform linear regression of node i on instrument i
        lm1 <- lm(data[, j] ~ instruments[, i])
        p_value_i_on_j <- summary(lm1)$coefficients[2, 4]  # p-value for instrument i on node i
        
        # Step 2.2: Perform linear regression of node j on instrument j
        lm2 <- lm(data[, i] ~ instruments[, j])
        p_value_j_on_i <- summary(lm2)$coefficients[2, 4]  # p-value for instrument j on node j
        
        # Check the conditions in Step 3
        if (p_value_i_on_j > level && p_value_j_on_i <= level) {
          is_root <- FALSE  # Node i cannot be a root
          break
        }
      }
    }
    
    # Step 4: Add i to the set of roots if it satisfies the conditions
    if (is_root) {
      root_set <- c(root_set, node_names[i])
    }
  }
  
  # Return the set of roots
  return(root_set)
}

root_nodes <- find_root_nodes(protein_network_sim, G_sim_fit)

################SORT FINDER
# Load necessary libraries
if (!requireNamespace("FNN", quietly = TRUE)) install.packages("FNN")
if (!requireNamespace("KernSmooth", quietly = TRUE)) install.packages("KernSmooth")
library(FNN)
library(KernSmooth)

# Function to perform linear regression and collect residuals
linear_regression_residuals <- function(data, instruments, pi) {
  residuals_list <- list()
  
  for (i in 1:ncol(instruments)) {
    if (!(colnames(instruments)[i] %in% pi)) {
      # Perform linear regression of the instrument on instruments in pi
      lm_fit <- lm(data[, i] ~ ., data = instruments[, pi, drop = FALSE])
      residuals_list[[colnames(data)[i]]] <- residuals(lm_fit)
    }
  }
  
  return(residuals_list)
}

# Function to compute pairwise linear regression between residuals
prune_unsorted <- function(residuals_list, U, level = 0.05) {
  pruned_U <- U  # Initialize pruned set as the unsorted set
  
  # Pairwise regression for residuals in U
  for (xi in U) {
    keep <- FALSE
    for (xj in setdiff(U, xi)) {
      lm1 <- lm(residuals_list[[xi]] ~ residuals_list[[xj]])
      lm2 <- lm(residuals_list[[xj]] ~ residuals_list[[xi]])
      
      # Get p-values for the slopes
      p_value_ij <- summary(lm1)$coefficients[2, 4]
      p_value_ji <- summary(lm2)$coefficients[2, 4]
      
      # Prune xi if the dependence condition is violated
      if (p_value_ij > level && p_value_ji <= level) {
        keep <- FALSE
        break
      } else {
        keep <- TRUE
      }
    }
    
    # Remove xi if it fails the dependence condition
    if (!keep) {
      pruned_U <- setdiff(pruned_U, xi)
    }
  }
  
  return(pruned_U)
}

# Function to compute mutual information (using correlation as a proxy)
mutual_information_proxy <- function(x, y) {
  cor_value <- cor(x, y, method = "pearson")
  return(-0.5 * log(1 - cor_value^2))  # Approximation based on correlation
}

# Function to compute the vertex with minimum dependence
find_min_dependence <- function(residuals_list, U, pi, data) {
  min_t_star <- Inf
  selected_node <- NULL
  
  for (node in U) {
    # Compute t*(node, pi)
    t_star_value <- sum(sapply(pi, function(xj) {
      mutual_information_proxy(residuals_list[[node]], data[,xj])
    }))
    
    if (t_star_value < min_t_star) {
      min_t_star <- t_star_value
      selected_node <- node
    }
  }
  
  return(selected_node)
}

# Sort Finder Algorithm
sort_finder <- function(data, instruments, roots, level = 0.05) {
  # Initialize π with roots, and unsorted set U as the remaining nodes
  pi <- roots
  U <- setdiff(colnames(data), pi)
  
  while (length(U) > 0) {
    # Stage 1: Prune U
    residuals_list <- linear_regression_residuals(data, instruments, pi)
    U_pruned <- prune_unsorted(residuals_list, U, level = 0.05)
    
    # Stage 2: Identify vertex with minimum dependence
    selected_node <- find_min_dependence(residuals_list, U_pruned, pi, data)
    if(is.null(selected_node)){
      break
    }
    # Add the selected node to π and remove it from U
    pi <- c(pi, selected_node)
    U <- setdiff(U, selected_node)
    
  }
  
  return(pi)
}



G_sim_fit1 = G_sim_fit
colnames(G_sim_fit1) = colnames(protein_network_sim)
sorted_nodes <- sort_finder(protein_network_sim, G_sim_fit1, root_nodes)
