#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// Gradient descent function for transfer learning
//' Compute beta joint
 //'
 //' @ NOTE: input data should be cleaned to standard format (SNP order should be same with reference SNP list)
 //' @param r0: summary statistics for base ancestry (set to 0 for ethnic-specific SNPs)
 //' @param n0: sample size for base ancestry 
 //' @param rk_list: list of summary statistics for auxillary ancestry (set to 0 for ethnic-specific SNPs)
 //' @param R0: correlation matrix of all SNPs for base ancestry (set to 0 for ethnic-specific SNPs)
 //' @param nk_list: list of sample size of auxillary ancestries
 //' @param Rk_list: correlation matrix of all SNPs for auxillary ancestry (set to 0 for ethnic-specific SNPs)
 //' @param indx: indicator of whether a SNP in a certain ethnic group
 //' @param M: number of ethnic groups
 //' @param alpha1-4: initialization parameters
 //' @return beta vector
// [[Rcpp::export]]
Rcpp::List gradient_descent_transfer_learning_rcpp_PRS(
    double n0,
    arma::vec r0,  // r0 is now a vector
    arma::mat R0,
    std::vector<double> nk_list,
    std::vector<arma::vec> rk_list,
    std::vector<arma::mat> Rk_list,
    double alpha1,
    double alpha2,
    double alpha3,
    double alpha4,
    double eta_l,
    double eta_m,
    int max_iter
) {
  
  int p = R0.n_rows;  // Assuming square matrices, p is the number of SNPs (rows in R0)
  
  // Initialize variables
  arma::vec u_l = alpha1 * arma::ones<arma::vec>(p);
  arma::vec v_l = alpha2 * arma::ones<arma::vec>(p);
  arma::vec h_m = alpha3 * arma::ones<arma::vec>(p);
  arma::vec g_m = alpha4 * arma::ones<arma::vec>(p);
  
  // Gradient descent for auxiliary data
  for (int l = 0; l <= max_iter; ++l) {
    arma::vec grad_u = arma::zeros<arma::vec>(p);
    arma::vec grad_v = arma::zeros<arma::vec>(p);
    
    for (size_t k = 0; k < nk_list.size(); ++k) {
      double nk = nk_list[k];
      arma::vec rk = rk_list[k];
      arma::mat Rk = Rk_list[k];
      
      arma::vec u_l_sq = arma::square(u_l);
      arma::vec v_l_sq = arma::square(v_l);
      arma::vec diff_sq = u_l_sq - v_l_sq;
      
      grad_u += (-4 * nk * rk % u_l + 4 * nk * (Rk * diff_sq) % u_l);
      grad_v += (4 * nk * rk % v_l - 4 * nk * (Rk * diff_sq) % v_l);
    }
    
    double total_n = std::accumulate(nk_list.begin(), nk_list.end(), 0.0) + n0;
    u_l -= (eta_l / total_n) * grad_u;
    v_l -= (eta_l / total_n) * grad_v;
    
    // Print the gradients for testing at each iteration (first few elements)
    if (l == 0 || l == max_iter / 2 || l == max_iter) {
      Rcpp::Rcout << "Iteration " << l << " - grad_u (first few): " << grad_u.head(5).t() << std::endl;
      Rcpp::Rcout << "Iteration " << l << " - grad_v (first few): " << grad_v.head(5).t() << std::endl;
    }
  }
  
  arma::vec hat_u = u_l;
  arma::vec hat_v = v_l;
  
  // Gradient descent for main data
  for (int m = 0; m <= max_iter; ++m) {
    arma::vec u_hat_sq = arma::square(hat_u);
    arma::vec v_hat_sq = arma::square(hat_v);
    arma::vec grad_h = arma::zeros<arma::vec>(p);
    arma::vec grad_g = arma::zeros<arma::vec>(p);
    
    arma::vec h_m_sq = arma::square(h_m);
    arma::vec g_m_sq = arma::square(g_m);
    arma::vec diff_sq = u_hat_sq - v_hat_sq + h_m_sq - g_m_sq;
    
    grad_h = (-4 * n0 * r0 % h_m + 4 * n0 * (R0 * diff_sq) % h_m);
    grad_g = (4 * n0 * r0 % g_m - 4 * n0 * (R0 * diff_sq) % g_m);
    
    h_m -= (eta_m / n0) * grad_h;
    g_m -= (eta_m / n0) * grad_g;
    
    // Print the gradients for testing at each iteration (first few elements)
    if (m == 0 || m == max_iter / 2 || m == max_iter) {
      Rcpp::Rcout << "Iteration " << m << " - grad_h (first few): " << grad_h.head(5).t() << std::endl;
      Rcpp::Rcout << "Iteration " << m << " - grad_g (first few): " << grad_g.head(5).t() << std::endl;
    }
  }
  
  arma::vec hat_h = h_m;
  arma::vec hat_g = g_m;
  
  // Define the output
  arma::vec hat_beta = arma::square(hat_u) - arma::square(hat_v) + arma::square(hat_h) - arma::square(hat_g);
  
  return Rcpp::List::create(
    Rcpp::Named("hat_u") = hat_u,
    Rcpp::Named("hat_v") = hat_v,
    Rcpp::Named("hat_h") = hat_h,
    Rcpp::Named("hat_g") = hat_g,
    Rcpp::Named("hat_beta") = hat_beta
  );
}
 
 
 // Function to apply gradient descent transfer learning across all blocks
 // [[Rcpp::export]]
 Rcpp::List gradient_descent_transfer_learning_all_blocks(
     Rcpp::List summ_list,
     Rcpp::List LD_list,
     int M,
     Rcpp::List indx,
     Rcpp::IntegerVector indx_block,
     double n0,
     std::vector<double> nk_list,
     double alpha1,
     double alpha2,
     double alpha3,
     double alpha4,
     double eta_l,
     double eta_m,
     int max_iter
 ) {
   int num_blocks = indx_block.size();
   Rcpp::List beta_results(num_blocks);
   
   // // Iterate over summ_list and replace NaNs with 0
   // for (size_t i = 0; i < summ_list.size(); ++i) {
   //   arma::mat summ = Rcpp::as<arma::mat>(summ_list[i]);
   //   summ.replace(arma::datum::nan, 0); // Replace NaN with 0
   //   summ_list[i] = summ; // Store back the modified matrix
   // }
   // 
   // // Iterate over LD_list and replace NaNs with 0
   // for (size_t i = 0; i < LD_list.size(); ++i) {
   //   arma::mat LD = Rcpp::as<arma::mat>(LD_list[i]);
   //   LD.replace(arma::datum::nan, 0); // Replace NaN with 0
   //   LD_list[i] = LD; // Store back the modified matrix
   // }
   
   for (int bl = 0; bl < num_blocks; ++bl) {
     if (indx_block[bl] == 0) {
       beta_results[bl] = Rcpp::List::create(Rcpp::Named("b") = R_NilValue);
       continue;
     }
     
     // Extract summary statistics and LD matrices
     arma::mat summ = Rcpp::as<arma::mat>(summ_list[bl]);
     Rcpp::List R = LD_list[bl];
     arma::mat indx_mat = Rcpp::as<arma::mat>(indx[bl]);
     
     // Prepare rk_list: all columns except the first from summ
     std::vector<arma::vec> rk_list;
     for (size_t i = 1; i < summ.n_cols; ++i) {
       rk_list.push_back(summ.col(i));
     }
     
     // Prepare Rk_list: all LD matrices except the first (target population)
     std::vector<arma::mat> Rk_list;
     for (size_t i = 1; i < R.size(); ++i) {
       Rk_list.push_back(Rcpp::as<arma::mat>(R[i]));
     }
     
     // Testing
     Rcpp::Rcout << "Checking inputs before calling gradient descent..." << std::endl;
     Rcpp::Rcout << "n0: " << n0 << ", total_n: " << std::accumulate(nk_list.begin(), nk_list.end(), 0.0) + n0 << std::endl;
     Rcpp::Rcout << "First few elements of r0: " << summ.col(0).head(10).t() << std::endl;
     Rcpp::Rcout << "First few elements of R0: " << Rcpp::as<arma::mat>(R[0]).submat(0,0,4,4) << std::endl;
     // Display first few elements of rk_list
     for (size_t i = 0; i < rk_list.size(); ++i) {
       Rcpp::Rcout << "First few elements of rk_list[" << i << "]: " << rk_list[i].head(10).t() << std::endl;
     }
     
     // Display first few elements of Rk_list
     for (size_t i = 0; i < Rk_list.size(); ++i) {
       Rcpp::Rcout << "First few elements of Rk_list[" << i << "] (top-left 5x5 block):\n" 
                   << Rk_list[i].submat(0, 0, std::min(4, static_cast<int>(Rk_list[i].n_rows) - 1),
       std::min(4, static_cast<int>(Rk_list[i].n_cols) - 1)) 
       << std::endl;
     }
     
     // Check dimensions of r0 (first column of summ)
     Rcpp::Rcout << "Dimension of r0 (summ.col(0)): " << summ.col(0).n_rows << " x " << summ.col(0).n_cols << std::endl;
     
     // Check dimensions of R0 (first LD matrix in the list)
     arma::mat R0 = Rcpp::as<arma::mat>(R[0]);
     Rcpp::Rcout << "Dimension of R0: " << R0.n_rows << " x " << R0.n_cols << std::endl;
     
     // Check dimensions of each rk_list element
     for (size_t i = 0; i < rk_list.size(); ++i) {
       Rcpp::Rcout << "Dimension of rk_list[" << i << "]: " << rk_list[i].n_rows << " x " << rk_list[i].n_cols << std::endl;
     }
     
     // Check dimensions of each Rk_list element
     for (size_t i = 0; i < Rk_list.size(); ++i) {
       Rcpp::Rcout << "Dimension of Rk_list[" << i << "]: " << Rk_list[i].n_rows << " x " << Rk_list[i].n_cols << std::endl;
     }
     
     // Call gradient descent function with correct inputs
     Rcpp::List beta_block = gradient_descent_transfer_learning_rcpp_PRS(
       n0, 
       summ.col(0), // r0 is the first column
       Rcpp::as<arma::mat>(R[0]), // R0 is the first LD matrix
       nk_list, 
       rk_list, 
       Rk_list,
       alpha1, alpha2, alpha3, alpha4, 
       eta_l, eta_m, max_iter
     );
     
     // Extract beta vector for the target population
     arma::vec beta_vec = as<arma::vec>(beta_block["hat_beta"]);
     Rcpp::Rcout << "First few elements of beta_block['hat_beta']: " << beta_vec.head(10).t() << std::endl;
     
     
     arma::vec indx_binary = arma::conv_to<arma::vec>::from(indx_mat.col(0) != 0);
     
     // Element-wise multiplication to apply indexing
     arma::vec beta_final = beta_vec % indx_binary;
     
     Rcpp::Rcout << "First few elements of beta_vec: " << beta_vec.head(10).t() << std::endl;
     
     // Store the final beta vector in the results
     beta_results[bl] = Rcpp::List::create(Rcpp::Named("b") = beta_final);
   }
   
   return beta_results;
 }
 