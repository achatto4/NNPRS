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
