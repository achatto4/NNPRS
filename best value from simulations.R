# Initialize an empty list named 'result'
result <- list()

# Loop through chromosome numbers 1 to 22
for (chr in 1:22) {
  file_path <- sprintf("/Users/anaghchattopadhyay/Documents/Nilanjan Project/MV-latent/chrom_res/chr%d_results.csv", chr)
  if (file.exists(file_path)) {
    result[[paste0("chr", chr)]] <- read.csv(file_path)
    print(paste("Loaded:", file_path))
  } else {
    print(paste("File not found:", file_path))
  }
}

# Initialize a dataframe to store the results
final_result <- data.frame(iter = integer(), eta = numeric(), alpha = numeric(), R2_sum = numeric())

# Loop through each chromosome's data
for (chr_data in result) {
  # Extract necessary columns
  temp_data <- chr_data[, c("iter", "eta", "alpha", "R2")]
  
  # Aggregate R2 by summing across matching iter, eta, and alpha
  final_result <- aggregate(R2 ~ iter + eta + alpha, data = rbind(final_result, temp_data), FUN = sum)
}

# Sort the final_result dataframe by R2 in descending order
final_result <- final_result[order(-final_result$R2), ]

##########
# Initialize an empty list to store the best rows
best_r2_rows <- list()

# Loop through each chromosome's data
for (i in 1:length(result)) {
  # Find the row with the maximum R2 value for this chromosome
  best_row <- result[[i]][which.max(result[[i]]$R2), ]
  
  # Store the best row in the list
  best_r2_rows[[i]] <- best_row
}

# Combine all the best rows into a matrix (or dataframe)
best_r2_matrix <- do.call(rbind, best_r2_rows)
##########

# Sum of the highest R2 values from all chromosomes
total_max_r2_sum <- sum(best_r2_matrix[,5])

###############

best_r2_matrix <- data.frame(
  chr = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22),
  iter = c(1, 2336, 546, 336, 127, 10000, 1, 2, 886, 30, 3, 1, 2336, 2, 78, 10000, 3, 3793, 4, 3, 6158, 4),
  eta = c(1e-02, 1e+00, 1e-04, 1e+00, 1e-04, 1e-02, 1e-01, 1e-03, 1e+00, 1e-04, 1e+00, 1e-04, 1e-02, 1e-01, 1e-04, 1e-01, 1e-01, 1e-01, 1e-01, 1e+00, 1e-01, 1e-02),
  alpha = c(0.0004882812, 0.0078125000, 0.1250000000, 0.1250000000, 0.0001220703, 0.5000000000, 0.5000000000, 0.0004882812, 0.0312500000, 0.0001220703, 0.1250000000, 0.0001220703, 0.5000000000, 0.5000000000, 0.0001220703, 0.0312500000, 0.5000000000, 0.0078125000, 0.5000000000, 0.1250000000, 0.1250000000, 0.5000000000)
)

# Initialize an empty list to store the best rows
best_r2_rows <- list()

# Loop through each chromosome's data
for (i in 1:length(result)) {
  # Find the row with the maximum R2 value for this chromosome
  best_row <- result[[i]][(result[[i]]$iter == best_r2_matrix[i,2] & result[[i]]$eta == best_r2_matrix[i,3] & result[[i]]$alpha == best_r2_matrix[i,4]),]
  
  # Store the best row in the list
  best_r2_rows[[i]] <- best_row
}

# Combine all the best rows into a matrix (or dataframe)
best_r2_matrix <- do.call(rbind, best_r2_rows)
best_r2_matrix

############rough

# Initialize an empty list to store the best rows
best_r2_rows <- list()

# Loop through each chromosome's data
for (i in 1:length(result)) {
  # Filter the data to include only rows where iter == 1 and eta == 1
  iter_eta_1_data <- result[[i]][result[[i]]$iter < 10, ]
  #  & result[[i]]$eta == 1 & result[[i]]$alpha == 1/2^13 
  # Find the row with the maximum R2 value for this chromosome where iter == 1 and eta == 1
  if (nrow(iter_eta_1_data) > 0) {
    best_row <- iter_eta_1_data[which.max(iter_eta_1_data$R2), ]
    
    # Store the best row in the list
    best_r2_rows[[i]] <- best_row
  }
}

# Combine all the best rows into a matrix (or dataframe)
best_r2_matrix_t <- do.call(rbind, best_r2_rows)
best_r2_matrix_t
